#!/usr/bin/env Rscript
# Find consensus peaks
FindConsensusPeaks <- function(gr, pwidth) {
  # 1. Define the target diameter
  # If pwidth is radius (150), target is 300
  target_width <- 2 * pwidth 
  
  # 2. Define the "footprint" of all peaks
  # Expand 1bp centers to full width to see where they touch
  gr_expanded <- resize(gr, width = target_width, fix = "center")
  
  # 3. Find "Islands" (Connected Components)
  # reduce merges anything that overlaps even by 1bp
  islands <- reduce(gr_expanded, min.gapwidth = 0L)
  
  # 4. Calculate Capacity
  # How many peaks fit in this island?
  # We use floor() to ensure we don't over-count and create overlaps. 
  # A 400bp island (target 300) gets 1 tile.
  island_capacities <- floor(width(islands) / target_width)
  
  # 5. Generate Tiles
  # tile() splits each island into 'n' equal-width adjacent chunks
  tiles <- unlist(tile(islands, n = island_capacities))
  
  # 6. Sanity Filter
  # Tiling is blind; it splits the space mathematically.
  # We must ensure each tile actually contains a real peak summit.
  # We overlap with the ORIGINAL 1bp centers (gr).
  tiles_with_peaks <- subsetByOverlaps(tiles, gr)
  
  return(tiles_with_peaks)
}

library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
}

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)

libs <- c("stf_2", "stf_3", "stf_4", "stf_5")
names(libs) <- libs

bedFiles_peak <- c(
  stf_2 = file.path("bam_new", "stf_2", "outs", "atac_peaks.bed"),
  stf_3 = file.path("bam_new", "stf_3", "outs", "atac_peaks.bed"),
  stf_4 = file.path("bam_new", "stf_4", "outs", "atac_peaks.bed"),
  stf_5 = file.path("bam_new", "stf_5", "outs", "atac_peaks.bed")
)

plist <- lapply(libs, function(x) {
  gr <- import(bedFiles_peak[[x]])
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- resize(gr, width = 1, fix = "center")
})

sgr <- Reduce(c, plist)
cpgr <- FindConsensusPeaks(sgr, 150)
export(cpgr, "consensus.bed")
