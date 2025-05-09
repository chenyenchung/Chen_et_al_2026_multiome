#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (length(args) != 3) {
  message(
    "\n",
    "Usage: \n",
    "--lib: The name of library as in the metadata\n",
    "--h5path: The path to cellranger-arc outputs. ",
    "Different libraries are expected to be subdirectories ",
    "under this path with subdirectory names specified by ",
    "--lib.\n",
    "--gtf: The GTF file used to generate the cellranger-arc ",
    "reference\n",
    "\n"
  )
  stop()
}

library(Seurat)
library(Signac)
library(Matrix)
library(rtracklayer)


if (interactive()) {
  h5path <- "../../../data/bam_new/stf_2/outs/raw_feature_bc_matrix.h5"
  lib <- "stf_2"
  gtfpath <- "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf"
} else {
  lib <- args$lib
  h5path <- file.path(args$h5, lib, "outs/raw_feature_bc_matrix.h5")
  gtfpath <- args$gtf
}

# Load matrices
matrices <- Read10X_h5(h5path)

# Get gene expression in bulk
bulk_exp <- rowSums(matrices[["Gene Expression"]])

# Get top genes
genes <- names(bulk_exp)[order(bulk_exp, decreasing = TRUE)]
top_genes <- genes[1:100]

# Get gene annotation
gtf <- import(gtfpath)
gtf <- subset(gtf, type == "gene" & gene_name %in% top_genes)

# Export as a BED
gtf$score <- 0
export(gtf, paste0(lib, "_mask.bed"))
