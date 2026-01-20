#!/usr/bin/env python
import argparse
import logging
import gzip
from collections import defaultdict
from pathlib import Path

logger = logging.getLogger('ins_depth')
logging.basicConfig(
  level=logging.DEBUG,
  format='[%(levelname)s]\t%(message)s'
)

parser = argparse.ArgumentParser(
  prog='Insertion Depth',
  description='Insertiong count per barcode'
)

parser.add_argument('--file', '-f', type=str, required=True)
parser.add_argument('--verbose', '-v', action='store_true')
parser.add_argument('--freq', '-q', type=int, default=1000000)
args = parser.parse_args()

input_path = Path(args.file)

depth = defaultdict(int)
iter = 0
with gzip.open(input_path, 'rt') as f:
  for line in f:
    iter += 1
    _seqname, _start, _rend, barcode, read_count = line.strip('\n').split('\t')
    depth[barcode] += int(read_count)

    if args.verbose and not iter % args.freq:
      logger.info(f'Processed {iter} entries...')


for bc, val in depth.items():
    print(f'{bc},{val}')
