#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (length(args) != 4) {
  message(
    "\n",
    "Usage:\n",
    "--lib: Library name as in the metadata\n",
    "--meta: Path to a CSV file containing 3 columns ",
    "(lib, line, domain) without a header\n",
    "--dgrp: Path to the driver line VCF file\n",
    "--out: Path to save the output\n",
    "\n"
  )
  stop()
}

if (interactive()) {
  this_lib <- "stf_2"
  meta_path <- "../../../data/line_used.txt"
  dgrp_vcf <- "/scratch/ycc520/thesis/work/d7/61af5a6227f0f33bf87304d8543dbc/dgrp_filtered.vcf.gz"
  outfile <- "../../../test.vcf.gz"
} else {
  this_lib <- args$lib
  meta_path <- args$meta
  dgrp_vcf <- args$dgrp
  outfile <- args$out
}

library(vcfR)

# Load metadata regarding which lines were used
# and their correspondence to the domain of origin
line_meta <- read.csv(meta_path, header = FALSE)
colnames(line_meta) <- c("lib", "line", "domain")

# Only keep the lines relevant to the current library
line_meta <- subset(line_meta, lib == this_lib)

# Load VCF files for DGRP variants
dgrp <- read.vcfR(dgrp_vcf)

# Generate output GT matrix
out_gt <- dgrp@gt[, c("FORMAT", line_meta$line)]
synth_sample <- paste(line_meta$line, line_meta$domain, sep = "@")
colnames(out_gt) <- c("FORMAT", synth_sample)

# We crossed it with a driver line, so F1 is expected to be a heterozygote
# if one of the parent (non-driver) carrying a variant.
for (icol in seq(2, ncol(out_gt))) {
  out_gt[, icol] <- ifelse(out_gt[, icol] == "1/1", "1/0", "0/0")
}

# Duplicate the vcfR object
sample_out <- dgrp

# Replace the GT matrix
sample_out@gt <- out_gt

write.vcf(sample_out, outfile)
