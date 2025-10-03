#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(R.utils))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (length(args) != 5) {
  message(
    "\n",
    "Usage:\n",
    "--lib: Library name as in the metadata\n",
    "--meta: Path to a CSV file containing 3 columns ",
    "(lib, line, domain) without a header\n",
    "--dgrp: Path to the DGRP VCF file\n",
    "--driver: Path to the driver line VCF file\n",
    "--out: Path to save the output\n",
    "\n"
  )
  stop()
}

if (interactive()) {
  this_lib <- "stf_2"
  meta_path <- "../../../data/line_used.txt"
  dgrp_vcf <- "../../../int/filtered_vcfs/dgrp_filtered.vcf.gz"
  driver_vcf <- "../../../int/filtered_vcfs/driver_filtered.vcf.gz"
  outfile <- "../../../test.vcf.gz"
} else {
  this_lib <- args$lib
  meta_path <- args$meta
  dgrp_vcf <- args$dgrp
  driver_vcf <- args$driver
  outfile <- args$out
}

library(vcfR)

# Load metadata regarding which lines were used
# and their correspondence to the domain of origin
line_meta <- read.csv(meta_path, header = FALSE)
colnames(line_meta) <- c("lib", "line", "domain")

# Only keep the lines relevant to the current library
line_meta <- subset(line_meta, lib == this_lib)

# Load VCF files for DGRP lines & driver-associated
# variants
dgrp <- read.vcfR(dgrp_vcf)
driver <- read.vcfR(driver_vcf)

# Subset only the lines used from the DGRP variant file
new_gt <- dgrp@gt[, c("FORMAT", line_meta$line)]

# Subset the GT field from the driver-associated variant file
dgt <- apply(driver@gt, 2, function(x) {
  # The field separator is :
  flist <- strsplit(x, ":")
  # The first field is GT
  gt <- vapply(flist, function(x) x[[1]], FUN.VALUE = character(1))
  return(gt)
})

# Copy the DGRP GT matrix for edit
out_gt <- new_gt

# Update the GT matrix by column
for (icol in seq(from = 2, to = ncol(out_gt))) {
  # Get the corresponding domain of origin from the metadata
  # Note that the first column of GT matrix is FORMAT, so there
  # is an offset of 1.
  domain <- line_meta$domain[icol - 1]

  # Get the GT value of the domain-specific driver
  thedriver <- dgt[, domain]

  # Get the GT value of the line
  theline <- new_gt[, icol]

  # Update the output GT matrix
  out_gt[, icol] <- mapply(
    # Take both GT info
    function(l, d) {
      # If the driver has a REF allele
      if (d %in% c("./.", ".|.", "0/0", "0|0")) {
        # Then the F1 genoypte is a HET or a homozygous REF
        out <- ifelse(l == "0/0", "0/0", "0/1")
      } else {
        # Otherwise the F1 genoypte is a HET or a homozygous ALT
        out <- ifelse(l == "0/0", "0/1", "1/1")
      }
      return(out)
    },
    l = theline,
    d = thedriver,
    SIMPLIFY = TRUE
  )
}

# Rename the samples with the domains of origin
synth_sample <- paste(line_meta$line, line_meta$domain, sep = "@")
colnames(out_gt) <- c("FORMAT", synth_sample)

# Duplicate the vcfR object
sample_out <- dgrp

# Replace the GT matrix
sample_out@gt <- out_gt

write.vcf(sample_out, outfile)
