#!/usr/bin/env Rscript
library(R.utils)
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$split),
  !is.null(args$lib),
  !is.null(args$file),
  file.exists(args$file)
)

renv::load(args$proot)

library(logger)
if (!is.null(args$debug)) {
  log_threshold(TRACE)
} else {
  log_threshold(INFO)
}

metadata <- read.csv(args$file, row.names = 1)
stopifnot("library" %in% colnames(metadata))
log_debug("Metadata columns: ", paste(colnames(metadata), collapse = ","))

if (grepl(",", args$split)) {
  split_cols <- strsplit(args$split, ",")[[1]]
} else {
  split_cols <- args$split
}

log_info("Split by: ", paste(split_cols, collapse = ", "))
stopifnot(all(split_cols %in% colnames(metadata)))

pre_output <- subset(metadata, library == args$lib)
stopifnot(nrow(pre_output) > 0)

out <- data.frame(
  barcode = sub("^stf_.#", "", row.names(pre_output)),
  cluster = do.call(function(...) paste(..., sep = "@"), pre_output[, split_cols])
)

write.csv(out, paste0(args$lib, ".csv"), row.names = FALSE, quote = FALSE)
