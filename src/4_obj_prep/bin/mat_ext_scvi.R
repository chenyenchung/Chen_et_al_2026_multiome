#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$obj <-"int/filtered_wnn.rds"
  args$proot <- "/scratch/ycc520/thesis"
  args$p <- "int/hvg_np/varp.rds"
  args$g <- "int/hvg_np/varg.rds"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$p),
  file.exists(args$p),
  !is.null(args$g),
  file.exists(args$g)
)

renv::load(args$proot)
library(Seurat)
library(Signac)
library(Matrix)

obj <- readRDS(args$obj)
varp <- readRDS(args$p)
varg <- readRDS(args$g)

p_tbl <- table(unlist(varp))
g_tbl <- table(unlist(varg))

g <- names(g_tbl[g_tbl > 2])
p <- names(p_tbl[p_tbl > 2])

gex <- LayerData(obj, layer = "counts", assay = "RNA")
gex <- as(gex, "RsparseMatrix")
gex <- gex[g, ]

writeLines(row.names(gex), "genes.txt")
writeLines(colnames(gex), "barcodes.txt")
writeMM(gex, "gex.mm")

peaks <- LayerData(obj, layer = "counts", assay = "ATAC")
peaks <- as(peaks, "RsparseMatrix")
peaks <- peaks[p, ]
writeLines(row.names(peaks), "peaks.txt")
writeMM(peaks, "peaks.mm")

write.csv(obj@meta.data, "metadata.csv")