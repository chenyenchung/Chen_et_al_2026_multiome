#!/usr/bin/env python
from scipy.io import mmread
import scanpy as sc
import pandas as pd
import anndata

pd.options.future.infer_string = False

# Load peak names and cell barcodes
peak_set = pd.read_table(
    'peaks.bed',
    sep='\t',
    header=None,
    index_col=False
)        
peak_set.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']

def make_peak_name(row):
    return f"{row['chr']}:{row['start']}-{row['end']}"
peak_set.index = peak_set.apply(make_peak_name, axis=1)

barcodes = []
with open('./barcodes.txt') as f:
    for line in f:
        barcodes.append(line.strip())
        

# Load metadata
metadata = pd.read_csv('./metadata.csv', index_col=0)

# Load ATAC peak counts
peak_counts = mmread('./peaks.mm').T.tocsr()

# Create AnnData object
peakdata = sc.AnnData(X=peak_counts)
peakdata.obs_names = barcodes
peakdata.var_names = peak_set.index.tolist()
peakdata.obs = metadata.loc[peakdata.obs_names]
peakdata.var = peak_set.loc[peakdata.var_names]


peakdata.write_h5ad('atac.h5ad')
