#!/usr/bin/env python
import argparse
import gzip
import logging
from pathlib import Path
import tempfile
import pyBigWig as pbw
from collections import defaultdict
import pandas as pd

logger = logging.getLogger("ins2bw")
logging.basicConfig(level=logging.DEBUG, format='[%(levelname)s]\t%(message)s')

parser = argparse.ArgumentParser(
    prog='Insertions to Bigwig',
    description='Convert 10x Genomics Fragment Files into BigWig Tracks'
)

parser.add_argument('--file', '-f', type=str, required=True)
parser.add_argument('--fai', '-g', type=str, required=True)
parser.add_argument('--metadata', '-m', type=str)
parser.add_argument('--split', '-i', type=str, required=True)
parser.add_argument('--chromosomes', '-c', type=str, default='X,Y,2R,2L,3R,3L,4')
parser.add_argument('--bin', '-b', type=int, default=20)
parser.add_argument('--scaling', '-s', type=int, default=1000000)
parser.add_argument('--depth', '-d', type=str, default='nCount_insertion')
parser.add_argument('--output', '-o', type=str, default='')

args = parser.parse_args()

input = Path(args.file)
if args.output != '':
    args.output = args.output + '_'

logger.debug(f'Commandline args: {args}')

# Get chromosome sizes
fai = Path(args.fai)
logger.info(f'Reading chromosome sizes from {args.fai}')
chrom_keep = args.chromosomes.split(',')
chrom_sizes = {}
with open(fai, 'rt') as f:
  for line in f:
    fields = line.strip('\n').split('\t')
    if fields[0] in chrom_keep:
      chrom_sizes[fields[0]] = int(fields[1])

logger.info(f'The following chromosomes are used: {", ".join(chrom_sizes)}')
logger.debug(chrom_sizes)

# Getting metadata if it exists
col_use = args.split.split(',')

metapath = Path(args.metadata)
meta = pd.read_csv(metapath, index_col=0)

cols_keep = []
for i in col_use:
    if i in meta.columns:
        cols_keep.append(i)
    else:
        logging.warning(f'{i} is not a column in metadata and will be ignored.')

logging.info(f"Generating bigwig files by {','.join(cols_keep)}")

if len(cols_keep) > 1:
    meta['split'] = meta[cols_keep[0]].str.cat(meta[cols_keep[1:]], sep='@')
else:
    meta['split'] = meta[cols_keep]

cluster_map = meta['split'].to_dict()

if args.depth in meta.columns:
    depth_map = meta[args.depth].to_dict()
else:
    logging.error(f'Depth column {args.depth} is not found in metadata.')
    raise ValueError
    

cluster_bins = defaultdict(lambda: defaultdict(int))
chroms_order = []

with gzip.open(input, 'rt') as f:
  for line in f:
    if line.startswith('#'): continue
    seqname, rstart, _rend, barcode, read_count = line.strip('\n').split('\t')

    if seqname not in chrom_sizes: continue
    if seqname not in chroms_order: chroms_order.append(seqname)

    rstart = int(rstart)
    read_count = int(read_count)

    # Determine the bin the read belongs to
    bin_start = (rstart // args.bin) * args.bin

    cluster = cluster_map[barcode]
    depth = depth_map[barcode]
    cluster_bins[cluster][(seqname, bin_start)] += read_count * args.scaling / depth

for cluster, bins in cluster_bins.items():
  # Sort by chromosome order, then position
  sorted_keys = sorted(
    bins.keys(),
    key=lambda k: (chroms_order.index(k[0]), k[1])
  )

  chroms = [k[0] for k in sorted_keys]
  starts = [k[1] for k in sorted_keys]
  ends = [min(k[1] + args.bin, chrom_sizes[k[0]]) for k in sorted_keys]
  values = [bins[k] for k in sorted_keys]

  with pbw.open(f'{args.output}{cluster}.bw', 'wb') as f:
    header = [(c, chrom_sizes[c]) for c in chroms_order]
    f.addHeader(header)
    f.addEntries(chroms, starts, ends=ends, values=values)
