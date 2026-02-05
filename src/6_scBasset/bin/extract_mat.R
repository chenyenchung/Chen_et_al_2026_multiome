#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$obj <- "int/objects/annotated.rds"
  args$proot <- "/scratch/ycc520/thesis"
  args$p <- "int/hvg_np/varp.rds"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$p),
  file.exists(args$p)
)

renv::load(args$proot)
library(Seurat)
library(Signac)
library(Matrix)
library(rtracklayer)
library(GenomicRanges)

obj <- readRDS(args$obj)
varp <- readRDS(args$p)

p_tbl <- table(unlist(varp))

p <- names(p_tbl[p_tbl > 2])
p <- p[!grepl("^Unmapped", p)]
p <- p[!grepl("^dmel", p)]
p <- p[!grepl("^rDNA", p)]

peaks <- LayerData(obj, layer = "counts", assay = "ATAC")
writeLines(colnames(peaks), "barcodes.txt")
peaks <- as(peaks, "RsparseMatrix")
peaks <- peaks[p, ]
export(GRanges(sub("-", ":", row.names(peaks))), "peaks.bed")
writeMM(peaks, "peaks.mm")

write.csv(obj@meta.data, "metadata.csv")
