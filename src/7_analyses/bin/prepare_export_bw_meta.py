#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


LIBRARY_RE = re.compile(r"^[A-Za-z0-9._-]+$")
CELL_RE = re.compile(r"^(?P<library>[^#]+)#(?P<barcode>.+)$")
MAX_LABEL_LENGTH = 63


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Normalize cluster-barcode CSV files for scbamop splitting."
    )
    parser.add_argument("--input", required=True, help="Cluster-barcode CSV file")
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Directory to write per-library scbamop metadata CSV files",
    )
    parser.add_argument(
        "--bam-root",
        required=True,
        help="Root directory containing <library>/outs/atac_possorted_bam.bam",
    )
    parser.add_argument(
        "--barcode-length",
        type=int,
        default=18,
        help="Expected cell barcode length after removing any library prefix",
    )
    return parser.parse_args()


def fail(message: str) -> None:
    raise SystemExit(f"ERROR: {message}")


def is_allowed_label_char(char: str) -> bool:
    return char.isascii() and (char.isalnum() or char in {"_", "-", " ", "."})


def sanitize_label(label: str) -> str:
    out = []
    for char in label:
        if char in {"/", "\\", "~"}:
            out.append("_")
        elif is_allowed_label_char(char):
            out.append(char)
        else:
            out.append("_")
    if out and out[0] == ".":
        out[0] = "_"
    sanitized = "".join(out)
    while ".." in sanitized:
        sanitized = sanitized.replace("..", "__")
    return sanitized


def get_label(row: dict[str, str], row_num: int) -> str:
    for key in ("cluster", "label"):
        if key in row and row[key] is not None:
            label = row[key].strip()
            if label:
                return label
    fail(f"row {row_num}: missing non-empty 'cluster' or 'label' column")


def split_cell(value: str, row_num: int) -> tuple[str, str]:
    match = CELL_RE.match(value)
    if match is None:
        fail(
            f"row {row_num}: expected a prefixed cell ID like 'stf_2#BARCODE' when "
            "library is not provided"
        )
    return match.group("library"), match.group("barcode")


def parse_row(row: dict[str, str], row_num: int) -> tuple[str, str, str]:
    library = row.get("library", "")
    library = library.strip() if library is not None else ""

    value = row.get("barcode")
    if value is None:
        value = row.get("cell")
    if value is None:
        fail(f"row {row_num}: missing 'barcode' or 'cell' column")
    value = value.strip()
    if not value:
        fail(f"row {row_num}: barcode or cell value is empty")

    if library:
        if "#" in value:
            parsed_library, barcode = split_cell(value, row_num)
            if parsed_library != library:
                fail(
                    f"row {row_num}: library column '{library}' does not match "
                    f"cell prefix '{parsed_library}'"
                )
        else:
            barcode = value
    else:
        library, barcode = split_cell(value, row_num)

    label = get_label(row, row_num)
    return library, barcode, label


def validate_library(library: str, row_num: int) -> None:
    if not library:
        fail(f"row {row_num}: library is empty")
    if not LIBRARY_RE.match(library):
        fail(
            f"row {row_num}: invalid library '{library}'. Use only letters, numbers, '.', '_', or '-'"
        )


def validate_barcode(barcode: str, row_num: int, barcode_length: int) -> None:
    if not barcode:
        fail(f"row {row_num}: barcode is empty")
    if "," in barcode:
        fail(f"row {row_num}: barcode must not contain commas")
    if len(barcode) != barcode_length:
        fail(
            f"row {row_num}: barcode length {len(barcode)} does not match expected "
            f"length {barcode_length}: {barcode}"
        )


def validate_label(label: str, row_num: int) -> None:
    if not label:
        fail(f"row {row_num}: cluster label is empty")
    if "," in label:
        fail(f"row {row_num}: cluster label must not contain commas")
    if len(label) > MAX_LABEL_LENGTH:
        fail(
            f"row {row_num}: cluster label is too long ({len(label)} chars; max {MAX_LABEL_LENGTH})"
        )


def validate_bams(bam_root: Path, libraries: list[str]) -> None:
    missing = []
    for library in libraries:
        bam_path = bam_root / library / "outs" / "atac_possorted_bam.bam"
        if not bam_path.exists():
            missing.append(str(bam_path))
    if missing:
        fail("missing ATAC BAM(s): " + ", ".join(missing))


def write_outputs(rows_by_library: dict[str, list[tuple[str, str]]], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for library, rows in rows_by_library.items():
        output_path = output_dir / f"{library}.csv"
        with output_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["barcode", "label"])
            writer.writerows(rows)


def main() -> None:
    args = parse_args()
    input_path = Path(args.input)
    output_dir = Path(args.output_dir)
    bam_root = Path(args.bam_root)

    if not input_path.exists():
        fail(f"input CSV not found: {input_path}")
    if not bam_root.exists():
        fail(f"BAM root directory not found: {bam_root}")

    rows_by_library: dict[str, list[tuple[str, str]]] = defaultdict(list)
    seen_barcodes: dict[str, set[str]] = defaultdict(set)
    sanitized_labels: dict[str, str] = {}

    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            fail("input CSV is missing a header row")

        required_any = {"cluster", "label"}
        if not required_any.intersection(reader.fieldnames):
            fail("input CSV must include either a 'cluster' or 'label' column")
        if not {"barcode", "cell"}.intersection(reader.fieldnames):
            fail("input CSV must include either a 'barcode' or 'cell' column")

        row_count = 0
        for row_num, row in enumerate(reader, start=2):
            if row is None:
                continue
            if not any(value and value.strip() for value in row.values()):
                continue

            library, barcode, label = parse_row(row, row_num)
            validate_library(library, row_num)
            validate_barcode(barcode, row_num, args.barcode_length)
            validate_label(label, row_num)

            sanitized = sanitize_label(label)
            existing = sanitized_labels.get(sanitized)
            if existing is None:
                sanitized_labels[sanitized] = label
            elif existing != label:
                fail(
                    f"row {row_num}: cluster label '{label}' collides with '{existing}' after "
                    f"scbamop sanitization to '{sanitized}'"
                )

            if barcode in seen_barcodes[library]:
                fail(
                    f"row {row_num}: duplicate barcode '{barcode}' detected in library '{library}'"
                )
            seen_barcodes[library].add(barcode)
            rows_by_library[library].append((barcode, label))
            row_count += 1

    if row_count == 0:
        fail("input CSV does not contain any usable records")

    libraries = list(rows_by_library)
    validate_bams(bam_root, libraries)
    write_outputs(rows_by_library, output_dir)

    print(
        f"Wrote {len(libraries)} metadata file(s) for {row_count} cell(s): "
        + ", ".join(libraries),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
