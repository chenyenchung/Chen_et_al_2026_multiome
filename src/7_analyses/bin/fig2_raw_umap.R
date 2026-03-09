#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "int/objects/annotated.rds"
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
library(ggplot2)
library(RColorBrewer)

obj <- readRDS(args$obj)

set.seed(1)
pal <- colorRampPalette(brewer.pal(n = 8, name = "Accent"))(43) |>
  sample()

p <- DimPlot(obj, group.by = "seurat_clusters") +
  scale_color_manual(
    values = pal
  ) +
  theme_void() +
  theme(
    plot.title = element_blank()
  )

if (!dir.exists("sup_fig/s2")) {dir.create("sup_fig/s2", recursive = TRUE)}
ggsave(
  filename = file.path("sup_fig", "s2", "seurat_cluster_umap.pdf"),
  plot = p, width = 6, height = 4
)
