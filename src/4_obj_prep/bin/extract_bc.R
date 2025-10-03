#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
if (interactive()) {
  args$obj <-"int/subset.rds"
  args$proot <- "/scratch/ycc520/thesis"
}
stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)
library(Seurat)
library(Signac)


obj_path <- args$obj
obj <- readRDS(obj_path)

cluster_bc <- split(row.names(obj[[]]), Idents(obj))
for (i in names(cluster_bc)) {
  barcodes <- strsplit(cluster_bc[[i]], "#")
  lib_barcodes <- split(
    sapply(barcodes, function(x) x[[2]]),
    sapply(barcodes, function(x) x[[1]])
  )
  for (j in names(lib_barcodes)) {
    writeLines(lib_barcodes[[j]], paste0(j, "_cluster_", i, ".txt"))
  }
}
