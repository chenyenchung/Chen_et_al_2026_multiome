#!/usr/bin/env python3
import argparse
import gzip
import json
import os
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import pysam
import tensorflow as tf

from scbasset.basenji_utils import dna_1hot
from scbasset.utils import make_model


NUCLEOTIDES = ["A", "C", "G", "T"]
SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def fail(message):
    raise SystemExit("ERROR: %s" % message)


def safe_token(value):
    token = SAFE_RE.sub("_", str(value).strip()).strip("_")
    if not token:
        token = "value"
    if token[0] == ".":
        token = "_" + token[1:]
    return token


def read_bed(path):
    rows = []
    with open(path) as handle:
        for line_num, raw in enumerate(handle, start=1):
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                parts = line.split("\t")
            else:
                parts = line.split(None, 3)
            if len(parts) < 3:
                fail("%s:%d has fewer than 3 BED columns" % (path, line_num))
            try:
                start = int(float(parts[1]))
                end = int(float(parts[2]))
            except ValueError:
                fail("%s:%d has non-numeric start/end coordinates" % (path, line_num))
            if end <= start:
                fail("%s:%d has end <= start" % (path, line_num))
            label = parts[3].strip() if len(parts) >= 4 and parts[3].strip() else "%s:%d-%d" % (parts[0], start, end)
            rows.append(
                {
                    "chrom": parts[0],
                    "start": start,
                    "end": end,
                    "label": label,
                    "line_num": line_num,
                }
            )
    if not rows:
        fail("BED file contains no usable records: %s" % path)
    return rows


def fetch_centered_sequence(fasta, chrom, start, end, seq_len):
    if chrom not in fasta.references:
        fail("chromosome '%s' is not present in FASTA %s" % (chrom, fasta.filename.decode() if isinstance(fasta.filename, bytes) else fasta.filename))

    mid = (start + end) // 2
    seq_start = mid - seq_len // 2
    seq_end = seq_start + seq_len
    fetch_start = max(0, seq_start)
    fetch_end = max(0, seq_end)

    prefix = "N" * max(0, -seq_start)
    seq = prefix + fasta.fetch(chrom, fetch_start, fetch_end).upper()
    if len(seq) < seq_len:
        seq += "N" * (seq_len - len(seq))
    if len(seq) > seq_len:
        seq = seq[:seq_len]
    return seq, seq_start, seq_end


def write_fasta_record(handle, name, seq, width=80):
    handle.write(">%s\n" % name)
    for i in range(0, len(seq), width):
        handle.write(seq[i : i + width] + "\n")


def command_prepare(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bed_rows = read_bed(args.bed)
    fasta = pysam.Fastafile(args.fasta)

    region_rows = []
    seen_ids = set()
    fasta_path = outdir / "regions.fa"
    bed_path = outdir / "regions_1344.bed"
    tsv_path = outdir / "regions.tsv"

    with fasta_path.open("w") as fa, bed_path.open("w") as bed_out:
        for idx, row in enumerate(bed_rows, start=1):
            region_id = "region_%04d" % idx
            if region_id in seen_ids:
                fail("duplicate generated region id: %s" % region_id)
            seen_ids.add(region_id)

            seq, seq_start, seq_end = fetch_centered_sequence(
                fasta, row["chrom"], row["start"], row["end"], args.seq_len
            )
            write_fasta_record(fa, region_id, seq)
            bed_out.write("%s\t%d\t%d\t%s\n" % (row["chrom"], max(0, seq_start), max(0, seq_end), region_id))
            region_rows.append(
                {
                    "region_id": region_id,
                    "label": row["label"],
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "seq_chrom": row["chrom"],
                    "seq_start": seq_start,
                    "seq_end": seq_end,
                    "seq_len": len(seq),
                    "bed_line": row["line_num"],
                }
            )

    fasta.close()
    pd.DataFrame(region_rows).to_csv(tsv_path, sep="\t", index=False)
    print("prepared %d region(s) in %s" % (len(region_rows), outdir), file=sys.stderr)


def read_fasta(path):
    seqs = {}
    current = None
    chunks = []
    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks).upper()
                current = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks).upper()
    if not seqs:
        fail("FASTA contains no records: %s" % path)
    return seqs


def parse_meme_text(path):
    motifs = {}
    current_id = None
    current_alt = None
    expected_width = None
    matrix = []

    def flush():
        if current_id is None:
            return
        if expected_width is not None and len(matrix) != expected_width:
            fail("motif %s expected width %d but parsed %d rows" % (current_id, expected_width, len(matrix)))
        if matrix:
            motifs[current_id] = {
                "motif_id": current_id,
                "motif_alt_id": current_alt,
                "pwm": np.asarray(matrix, dtype=np.float32),
            }

    with open(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("MOTIF "):
                flush()
                parts = line.split()
                current_id = parts[1]
                current_alt = parts[2] if len(parts) > 2 else ""
                expected_width = None
                matrix = []
                continue
            if line.startswith("letter-probability matrix:"):
                match = re.search(r"\bw=\s*(\d+)", line)
                expected_width = int(match.group(1)) if match else None
                matrix = []
                continue
            if current_id is not None and expected_width is not None and len(matrix) < expected_width:
                vals = line.split()
                if len(vals) >= 4:
                    try:
                        matrix.append([float(vals[0]), float(vals[1]), float(vals[2]), float(vals[3])])
                    except ValueError:
                        pass
    flush()
    if not motifs:
        fail("no motifs parsed from MEME file: %s" % path)
    return motifs


def read_fimo(path):
    if not Path(path).exists():
        fail("FIMO TSV not found: %s" % path)
    try:
        hits = pd.read_csv(path, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        hits = pd.DataFrame()
    return hits


def bh_qvalues(pvalues):
    pvalues = np.asarray(pvalues, dtype=np.float64)
    qvalues = np.ones_like(pvalues)
    if pvalues.size == 0:
        return qvalues
    order = np.argsort(pvalues)
    ranked = pvalues[order] * float(pvalues.size) / np.arange(1, pvalues.size + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvalues[order] = np.minimum(ranked, 1.0)
    return qvalues


def make_sigmoid_predictor(model):
    latent_model = tf.keras.Model(inputs=model.layers[0].input, outputs=model.layers[-4].output)
    weights = model.layers[-3].get_weights()[0]

    def predict(x):
        x_tf = tf.convert_to_tensor(x, dtype=tf.float32)
        n_sequences = int(x_tf.shape[0])
        latent = latent_model(x_tf, training=False).numpy()
        latent = np.asarray(latent)
        try:
            latent = latent.reshape((n_sequences, -1))
        except ValueError:
            fail(
                "latent model output shape %s cannot be reshaped to %d sequence(s)"
                % (latent.shape, n_sequences)
            )
        if latent.shape[1] != weights.shape[0]:
            fail(
                "latent model output has %d feature(s), expected %d for classifier weights"
                % (latent.shape[1], weights.shape[0])
            )
        logits = np.dot(latent, weights)
        pred = 1.0 / (1.0 + np.exp(-logits))
        if pred.ndim != 2:
            fail("prediction output has shape %s, expected a 2D matrix" % (pred.shape,))
        return pred

    return predict


def make_groups(obs, pair):
    columns = [c.strip() for c in pair.split(",") if c.strip()]
    if len(columns) != 2:
        fail("--groupby-pair must contain exactly two comma-separated columns, got: %s" % pair)
    for col in columns:
        if col not in obs.columns:
            fail("metadata column '%s' not found in ad.obs" % col)

    valid_mask = ~obs.loc[:, columns].isna().any(axis=1)
    n_valid = int(valid_mask.sum())
    n_excluded = int((~valid_mask).sum())
    if n_valid < 2:
        fail("grouping '%s' has fewer than two cells after excluding missing metadata" % pair)
    if n_excluded:
        print(
            "excluding %d cell(s) with missing values for grouping '%s'" % (n_excluded, pair),
            file=sys.stderr,
        )

    values = obs.loc[valid_mask, columns].astype(str)
    labels = values[columns[0]].str.cat(values[columns[1]], sep=" x ")
    unique_labels = sorted(labels.unique())
    if len(unique_labels) < 2:
        fail("grouping '%s' produced fewer than two groups" % pair)

    label_to_id = {label: i for i, label in enumerate(unique_labels)}
    group_index = labels.map(label_to_id).values.astype(np.int64)
    sizes = np.bincount(group_index, minlength=len(unique_labels)).astype(np.float32)
    if np.any(sizes == 0):
        fail("grouping '%s' produced an empty group" % pair)
    if np.any(sizes == float(n_valid)):
        fail("grouping '%s' has a group containing all cells, so group-vs-rest is undefined" % pair)

    indicator = np.zeros((obs.shape[0], len(unique_labels)), dtype=np.float32)
    valid_positions = np.where(valid_mask.values)[0]
    indicator[valid_positions, group_index] = 1.0
    valid_vector = valid_mask.values.astype(np.float32)

    group_rows = []
    used_ids = set()
    for label in unique_labels:
        parts = label.split(" x ", 1)
        group_id = safe_token("%s__%s" % (parts[0], parts[1]))
        base_id = group_id
        counter = 2
        while group_id in used_ids:
            group_id = "%s_%d" % (base_id, counter)
            counter += 1
        used_ids.add(group_id)
        idx = label_to_id[label]
        group_rows.append(
            {
                "group_id": group_id,
                "group_label": label,
                columns[0]: parts[0],
                columns[1]: parts[1],
                "n_cells": int(sizes[idx]),
                "n_rest": int(n_valid - sizes[idx]),
            }
        )

    return columns, labels, indicator, valid_vector, sizes, pd.DataFrame(group_rows)


def aggregate_predictions(pred, indicator, valid_vector, sizes):
    if pred.ndim != 2:
        fail("prediction output has shape %s, expected (batch, cells)" % (pred.shape,))
    if pred.shape[1] != indicator.shape[0]:
        fail(
            "prediction output has %d cell(s), expected %d from AnnData"
            % (pred.shape[1], indicator.shape[0])
        )
    sums = np.dot(pred, indicator)
    total_valid = float(valid_vector.sum())
    total = np.dot(pred, valid_vector[:, np.newaxis])
    group_mean = sums / sizes[np.newaxis, :]
    rest_mean = (total - sums) / (total_valid - sizes)[np.newaxis, :]
    return group_mean, rest_mean


def compute_region_ism(seq, model_predict, indicator, valid_vector, sizes, seq_len, batch_size):
    seq_ref = dna_1hot(seq, seq_len=seq_len, n_uniform=True).astype(np.float32)
    ref_pred = model_predict(seq_ref[np.newaxis, :, :])
    ref_group, ref_rest = aggregate_predictions(ref_pred, indicator, valid_vector, sizes)

    n_groups = indicator.shape[1]
    group_ism = np.zeros((n_groups, seq_len, 4), dtype=np.float32)
    rest_ism = np.zeros((n_groups, seq_len, 4), dtype=np.float32)
    mutation_ids = np.arange(seq_len * 4)

    for start in range(0, mutation_ids.size, batch_size):
        ids = mutation_ids[start : start + batch_size]
        batch = np.repeat(seq_ref[np.newaxis, :, :], len(ids), axis=0)
        positions = ids // 4
        nts = ids % 4
        batch[np.arange(len(ids)), positions, :] = 0.0
        batch[np.arange(len(ids)), positions, nts] = 1.0

        pred = model_predict(batch)
        group_mean, rest_mean = aggregate_predictions(pred, indicator, valid_vector, sizes)
        group_ism[:, positions, nts] = (group_mean - ref_group).T
        rest_ism[:, positions, nts] = (rest_mean - ref_rest).T

    group_ism -= group_ism.mean(axis=2, keepdims=True)
    rest_ism -= rest_ism.mean(axis=2, keepdims=True)
    specific_ism = group_ism - rest_ism
    group_ref = (group_ism * seq_ref[np.newaxis, :, :]).sum(axis=2)
    rest_ref = (rest_ism * seq_ref[np.newaxis, :, :]).sum(axis=2)
    specific_ref = group_ref - rest_ref
    return group_ism, rest_ism, specific_ism, group_ref, rest_ref, specific_ref


def dot_pwm(pwm, ism_matrix):
    width = min(pwm.shape[0], ism_matrix.shape[0])
    if width <= 0:
        return np.nan
    return float(np.dot(pwm[:width, :].reshape(-1), ism_matrix[:width, :].reshape(-1)))


def sliding_pwm_scores(pwm, track):
    width = pwm.shape[0]
    if width <= 0 or track.shape[0] < width:
        return np.asarray([], dtype=np.float64)
    scores = np.zeros(track.shape[0] - width + 1, dtype=np.float64)
    for nt_idx in range(4):
        scores += np.correlate(track[:, nt_idx], pwm[:, nt_idx], mode="valid")
    return scores


def write_table(df, path):
    path = Path(path)
    if path.suffix == ".gz":
        with gzip.open(path, "wt") as handle:
            df.to_csv(handle, sep="\t", index=False)
    else:
        df.to_csv(path, sep="\t", index=False)


def command_score(args):
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    tf.get_logger().setLevel("ERROR")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    regions = pd.read_csv(args.regions_tsv, sep="\t")
    required_region_cols = {"region_id", "label", "chrom", "start", "end", "seq_chrom", "seq_start", "seq_end", "seq_len"}
    missing = required_region_cols.difference(regions.columns)
    if missing:
        fail("regions TSV is missing columns: %s" % ", ".join(sorted(missing)))

    seqs = read_fasta(args.fasta)
    for region_id in regions["region_id"]:
        if region_id not in seqs:
            fail("region %s is missing from FASTA %s" % (region_id, args.fasta))
        if len(seqs[region_id]) != args.seq_len:
            fail("region %s has FASTA length %d, expected %d" % (region_id, len(seqs[region_id]), args.seq_len))

    import anndata

    ad = anndata.read_h5ad(args.anndata)
    columns, group_labels, indicator, valid_vector, sizes, groups = make_groups(ad.obs, args.groupby_pair)
    groups.to_csv(outdir / "groups.tsv", sep="\t", index=False)
    regions.to_csv(outdir / "regions.tsv", sep="\t", index=False)

    model = make_model(args.bottleneck, ad.shape[0], seq_len=args.seq_len, show_summary=False)
    model.load_weights(args.model)
    predictor = make_sigmoid_predictor(model)

    n_regions = regions.shape[0]
    n_groups = groups.shape[0]
    h5_path = outdir / "group_ism.h5"
    str_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(h5_path, "w") as h5:
        h5.attrs["seq_len"] = args.seq_len
        h5.attrs["groupby_pair"] = args.groupby_pair
        h5.attrs["groupby_columns"] = json.dumps(columns)
        h5.create_dataset("region_ids", data=regions["region_id"].astype(str).values, dtype=str_dtype)
        h5.create_dataset("region_labels", data=regions["label"].astype(str).values, dtype=str_dtype)
        h5.create_dataset("group_ids", data=groups["group_id"].astype(str).values, dtype=str_dtype)
        h5.create_dataset("group_labels", data=groups["group_label"].astype(str).values, dtype=str_dtype)
        h5.create_dataset("nucleotides", data=np.asarray(NUCLEOTIDES, dtype=object), dtype=str_dtype)
        kwargs4 = dict(shape=(n_regions, n_groups, args.seq_len, 4), dtype="float32", compression="gzip")
        kwargs3 = dict(shape=(n_regions, n_groups, args.seq_len), dtype="float32", compression="gzip")
        ds_group = h5.create_dataset("group_ism", **kwargs4)
        ds_rest = h5.create_dataset("rest_ism", **kwargs4)
        ds_specific = h5.create_dataset("specific_ism", **kwargs4)
        ds_group_ref = h5.create_dataset("group_ref", **kwargs3)
        ds_rest_ref = h5.create_dataset("rest_ref", **kwargs3)
        ds_specific_ref = h5.create_dataset("specific_ref", **kwargs3)

        region_arrays = {}
        for r_idx, region in regions.iterrows():
            region_id = region["region_id"]
            print("computing ISM for %s (%d/%d)" % (region_id, r_idx + 1, n_regions), file=sys.stderr)
            arrays = compute_region_ism(
                seqs[region_id],
                predictor,
                indicator,
                valid_vector,
                sizes,
                args.seq_len,
                args.mutation_batch_size,
            )
            group_ism, rest_ism, specific_ism, group_ref, rest_ref, specific_ref = arrays
            ds_group[r_idx, :, :, :] = group_ism
            ds_rest[r_idx, :, :, :] = rest_ism
            ds_specific[r_idx, :, :, :] = specific_ism
            ds_group_ref[r_idx, :, :] = group_ref
            ds_rest_ref[r_idx, :, :] = rest_ref
            ds_specific_ref[r_idx, :, :] = specific_ref
            region_arrays[region_id] = {
                "group_ism": group_ism,
                "rest_ism": rest_ism,
                "specific_ism": specific_ism,
                "specific_ref": specific_ref,
            }

    motifs = parse_meme_text(args.meme)
    fimo = read_fimo(args.fimo_tsv)
    motif_rows = build_motif_hits(fimo, motifs, regions, groups, region_arrays)
    motif_hits = pd.DataFrame(motif_rows)
    if not motif_hits.empty:
        motif_hits["empirical_qvalue"] = bh_qvalues(motif_hits["empirical_pvalue"].values)
        motif_hits["abs_specificity_score"] = motif_hits["specificity_score"].abs()
        motif_hits = motif_hits.sort_values(
            ["group_id", "region_id", "abs_specificity_score"],
            ascending=[True, True, False],
        )
        motif_hits["rank_in_group"] = motif_hits.groupby("group_id")["abs_specificity_score"].rank(
            method="first", ascending=False
        ).astype(int)
    write_table(motif_hits, outdir / "motif_hits.tsv.gz")

    top_windows = pd.DataFrame(build_top_windows(regions, groups, region_arrays, seqs, args.window_size, args.top_windows))
    write_table(top_windows, outdir / "top_windows.tsv.gz")

    summary = {
        "anndata": str(args.anndata),
        "model": str(args.model),
        "regions_tsv": str(args.regions_tsv),
        "fasta": str(args.fasta),
        "fimo_tsv": str(args.fimo_tsv),
        "meme": str(args.meme),
        "groupby_pair": args.groupby_pair,
        "groupby_columns": columns,
        "n_cells": int(ad.shape[0]),
        "n_grouped_cells": int(valid_vector.sum()),
        "n_excluded_cells": int(ad.shape[0] - valid_vector.sum()),
        "n_regions": int(n_regions),
        "n_groups": int(n_groups),
        "n_fimo_rows": int(fimo.shape[0]) if not fimo.empty else 0,
        "n_motif_hit_group_rows": int(motif_hits.shape[0]),
        "n_top_window_rows": int(top_windows.shape[0]),
        "seq_len": int(args.seq_len),
        "mutation_batch_size": int(args.mutation_batch_size),
    }
    with (outdir / "run_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


def build_motif_hits(fimo, motifs, regions, groups, region_arrays):
    if fimo.empty:
        return []
    required_cols = {"motif_id", "sequence_name", "start", "stop", "strand"}
    missing = required_cols.difference(fimo.columns)
    if missing:
        fail("FIMO TSV is missing columns: %s" % ", ".join(sorted(missing)))

    region_meta = regions.set_index("region_id").to_dict(orient="index")
    rows = []
    background_cache = {}
    for _, hit in fimo.iterrows():
        region_id = str(hit["sequence_name"])
        if region_id not in region_arrays:
            continue
        motif_id = str(hit["motif_id"])
        if motif_id not in motifs:
            fail("FIMO motif_id '%s' was not found in MEME file" % motif_id)

        pwm = motifs[motif_id]["pwm"]
        strand = str(hit["strand"])
        if strand == "-":
            pwm = pwm[::-1, ::-1]
        start0 = int(hit["start"]) - 1
        stop0 = int(hit["stop"])
        width = min(pwm.shape[0], stop0 - start0, region_arrays[region_id]["specific_ism"].shape[1] - start0)
        if start0 < 0 or width <= 0:
            continue
        pwm = pwm[:width, :]

        region = region_meta[region_id]
        motif_genomic_start = int(region["seq_start"]) + start0
        motif_genomic_end = int(region["seq_start"]) + start0 + width

        for g_idx, group in groups.iterrows():
            group_matrix = region_arrays[region_id]["group_ism"][g_idx, start0 : start0 + width, :]
            rest_matrix = region_arrays[region_id]["rest_ism"][g_idx, start0 : start0 + width, :]
            specific_track = region_arrays[region_id]["specific_ism"][g_idx, :, :]
            specific_matrix = specific_track[start0 : start0 + width, :]

            group_score = dot_pwm(pwm, group_matrix)
            rest_score = dot_pwm(pwm, rest_matrix)
            cache_key = (region_id, int(g_idx), motif_id, strand, width)
            bg_scores = background_cache.get(cache_key)
            if bg_scores is None:
                bg_scores = sliding_pwm_scores(pwm, specific_track)
                background_cache[cache_key] = bg_scores
            if start0 < bg_scores.size:
                specificity_score = float(bg_scores[start0])
            else:
                specificity_score = dot_pwm(pwm, specific_matrix)
            pvalue = float((np.sum(np.abs(bg_scores) >= abs(specificity_score)) + 1.0) / (bg_scores.size + 1.0)) if bg_scores.size else 1.0

            out = hit.to_dict()
            out.update(
                {
                    "region_id": region_id,
                    "region_label": region["label"],
                    "region_chrom": region["chrom"],
                    "region_start": int(region["start"]),
                    "region_end": int(region["end"]),
                    "motif_genomic_start": motif_genomic_start,
                    "motif_genomic_end": motif_genomic_end,
                    "motif_width": width,
                    "group_id": group["group_id"],
                    "group_label": group["group_label"],
                    "n_cells": int(group["n_cells"]),
                    "n_rest": int(group["n_rest"]),
                    "group_score": group_score,
                    "rest_score": rest_score,
                    "specificity_score": specificity_score,
                    "score_per_bp": specificity_score / float(width),
                    "direction": "group_higher" if specificity_score >= 0 else "rest_higher",
                    "empirical_pvalue": pvalue,
                }
            )
            rows.append(out)
    return rows


def build_top_windows(regions, groups, region_arrays, seqs, window_size, top_n):
    if top_n <= 0:
        return []
    rows = []
    region_meta = regions.set_index("region_id").to_dict(orient="index")
    for region_id, arrays in region_arrays.items():
        region = region_meta[region_id]
        seq = seqs[region_id]
        specific_ref = arrays["specific_ref"]
        if specific_ref.shape[1] < window_size:
            continue
        for g_idx, group in groups.iterrows():
            track = specific_ref[g_idx, :]
            scores = np.array(
                [track[start : start + window_size].mean() for start in range(0, track.shape[0] - window_size + 1)],
                dtype=np.float64,
            )
            order = np.argsort(np.abs(scores))[::-1][:top_n]
            for rank, start0 in enumerate(order, start=1):
                end0 = int(start0) + window_size
                rows.append(
                    {
                        "region_id": region_id,
                        "region_label": region["label"],
                        "group_id": group["group_id"],
                        "group_label": group["group_label"],
                        "rank": rank,
                        "window_start_0based": int(start0),
                        "window_end_0based": end0,
                        "genomic_chrom": region["seq_chrom"],
                        "genomic_start": int(region["seq_start"]) + int(start0),
                        "genomic_end": int(region["seq_start"]) + end0,
                        "sequence": seq[int(start0) : end0],
                        "specificity_ref_mean": float(scores[int(start0)]),
                        "abs_specificity_ref_mean": float(abs(scores[int(start0)])),
                    }
                )
    return rows


def make_parser():
    parser = argparse.ArgumentParser(description="Run scBasset ISM and motif scoring for BED regions.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prep = subparsers.add_parser("prepare", help="prepare centered FASTA/BED records from input BED")
    prep.add_argument("--bed", required=True)
    prep.add_argument("--fasta", required=True)
    prep.add_argument("--outdir", required=True)
    prep.add_argument("--seq-len", type=int, default=1344)
    prep.set_defaults(func=command_prepare)

    score = subparsers.add_parser("score", help="score prepared regions with scBasset ISM")
    score.add_argument("--regions-tsv", required=True)
    score.add_argument("--fasta", required=True)
    score.add_argument("--fimo-tsv", required=True)
    score.add_argument("--meme", required=True)
    score.add_argument("--anndata", required=True)
    score.add_argument("--model", required=True)
    score.add_argument("--groupby-pair", required=True)
    score.add_argument("--outdir", required=True)
    score.add_argument("--seq-len", type=int, default=1344)
    score.add_argument("--bottleneck", type=int, default=32)
    score.add_argument("--mutation-batch-size", type=int, default=256)
    score.add_argument("--top-windows", type=int, default=100)
    score.add_argument("--window-size", type=int, default=10)
    score.set_defaults(func=command_score)
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
