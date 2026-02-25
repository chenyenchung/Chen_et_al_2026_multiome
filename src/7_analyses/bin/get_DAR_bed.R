#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$sdar <- "result/de_spatial/ATAC.csv"
}

stopifnot(
  !is.null(args$sdar),
  file.exists(args$sdar),
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)
library(rtracklayer)

dar_tbl <- read.csv(args$sdar)
sig_dar <- subset(dar_tbl, FDR < 0.05)

peak_list <- split(
  sig_dar$feature,
  paste(sig_dar$contrast, vapply(
    sig_dar$logFC > 0, function(x) ifelse(x, "UP", "DOWN"),
    FUN.VALUE = character(1)
  ))
)
peak_list <- lapply(peak_list, unique)

gr_list <- lapply(peak_list, function(peaks) {
  peak_ele <- strsplit(peaks, "-")
  peak_f <- vapply(peak_ele, function(x) {
    x_len <- length(x)
    coords <- paste(x[seq(x_len - 1, x_len)], collapse = "-")
    sname <- paste(x[seq(1, x_len - 2)], collapse = "_")
    out <- paste(sname, coords, sep = ":")
    return(out)
  }, FUN.VALUE = character(1))
  gr <- GRanges(peak_f)
  gr <- sortSeqlevels(gr)
  gr <- sort(gr)
  return(gr)
})

for (bed in names(gr_list)) {
  out_fn <- paste0(gsub(" ", "-", bed), ".bed")
  export(gr_list[[bed]], out_fn)
}
