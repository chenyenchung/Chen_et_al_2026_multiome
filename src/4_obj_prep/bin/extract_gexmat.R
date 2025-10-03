#!/usr/bin/env Rscript
library(R.utils)
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "int/base_obj.rds"
  args$mkrds <- "/scratch/ycc520/thesis/data/NN_asset/MarkersNNP15.rds"
}

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$mkrds),
  file.exists(args$mkrds)
)


renv::load(args$proot)
library(Seurat)
library(Signac)

obj <- readRDS(args$obj)

if (!"data" %in% Layers(obj, assay = "RNA")) {
  obj <- NormalizeData(obj, assay = "RNA")
}
outmat <- LayerData(obj, layer = "data", assay = "RNA")
outmat <- as(outmat, "RsparseMatrix")

# Get P15 markers
mks <- readRDS(args$mkrds)

# Pad markers if necessary
if (!all(mks %in% row.names(outmat))) {
  message(
    "Padding the expression matrix with rows of 0 for undetected markers..."
  )
  missing_mks <- setdiff(mks, row.names(outmat))
  pad <- matrix(0, nrow = length(missking_mks), ncol = ncol(outmat))
  colnames(pad) <- colnames(outmat)
  row.names(pad) <- missing_mks
  pad <- as(pad, "RsparseMatrix")
  outmat <- Matrix::rbind2(outmat, pad)
}

outmat <- outmat[mks, ]
outmat <- as.matrix(outmat)

saveRDS(outmat, "gex_mat.rds")
