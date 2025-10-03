#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

print(args)
if (length(args) != 2) {
  message(
    "\n",
    "Usage: \n",
    "--lib: The name of library as in the metadata\n",
    "--h5path: The path to cellranger-arc outputs. ",
    "Different libraries are expected to be subdirectories ",
    "under this path with subdirectory names specified by ",
    "--lib.\n",
    "\n"
  )
  stop()
}

library(Seurat)
library(Signac)
library(Matrix)


if (interactive()) {
  h5path <- "../../data/bam_new/stf_2/outs/raw_feature_bc_matrix.h5"
  lib <- "stf_2"
} else {
  lib <- args$lib
  h5path <- file.path(args$h5, lib, "outs/raw_feature_bc_matrix.h5")
}

# Load matrices
matrices <- Read10X_h5(h5path)

# Gene expression cutoff
g2keep <- colSums(matrices[["Gene Expression"]] > 0) > 200L

# ATAC cutoff
p2keep <- colSums(matrices[["Peaks"]]) > 300L

# Accepted barcodes
abc <- colnames(matrices[["Gene Expression"]])[g2keep & p2keep]

# Save barcode list for demuxlet
writeLines(abc, paste(lib, "txt", sep = "."))
