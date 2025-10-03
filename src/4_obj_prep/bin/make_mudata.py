#!/usr/bin/env python
import os
from scipy.io import mmread
import scanpy as sc
import pandas as pd
import muon as mu
import anndata as adata

# Load RNA counts
rna_counts = mmread('./gex.mm').T.tocsr()

# Load gene and cell barcodes
genes = []
with open('./genes.txt') as f:
    for line in f:
        genes.append(line.strip())
        
barcodes = []
with open('./barcodes.txt') as f:
    for line in f:
        barcodes.append(line.strip())
        

# Load metadata
metadata = pd.read_csv('./metadata.csv', index_col=0)

# Create AnnData object
adata = sc.AnnData(X=rna_counts)
adata.var_names = genes
adata.obs_names = barcodes

# Load ATAC peak counts
peak_counts = mmread('./peaks.mm').T.tocsr()

# Load peak names
peaks = []
with open('./peaks.txt') as f:
    for line in f:
        peaks.append(line.strip())

# Create AnnData object
peakdata = sc.AnnData(X=peak_counts)
peakdata.var_names = peaks
peakdata.obs_names = barcodes

mdata = mu.MuData({'rna': adata, 'atac': peakdata})
mdata.obs = metadata.loc[mdata.obs_names]
print(mdata)

mdata.write_h5mu('./obj.h5mu')
