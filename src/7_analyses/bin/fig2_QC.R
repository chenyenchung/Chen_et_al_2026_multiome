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
library(patchwork)
library(data.table)

obj <- readRDS(args$obj)
metadata <- obj@meta.data |>
  as.data.table()

lib_size <- metadata[, .N, by = "library"]

p_lib <- lib_size |>
  ggplot(aes(x = library, y = N, fill = library)) +
  geom_col() +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  theme_classic() +
  labs(
    x = "Library",
    y = "# of Cells"
  )

slib <- metadata[, .N, by = c("library", "spatial_origin")][
  , lib_size := sum(N), by = "library"
][, prop := N / lib_size]

slib$spatial_origin <- factor(
  slib$spatial_origin,
  levels = c("pxb", "optix", "dpp"),
  labels = c("Vsx", "Optix", "Bi")
)

p_spa <- slib |>
  ggplot(aes(x = library, y = prop, fill = spatial_origin)) +
  geom_col() +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(
    values = c(
      "Vsx" = "#1B9E77",
      "Optix" = "#7570B3",
      "Bi" = "#D95F02"
    )
  ) +
  theme_classic() +
  labs(
    x = "Library",
    y = "% of Library",
    fill = "Spatial Origin"
  )

tlib <- metadata[, .N, by = c("library", "temporal_identity")][
  , lib_size := sum(N), by = "library"
][, prop := N / lib_size]
tlib$temporal_identity <- factor(
  tlib$temporal_identity,
  levels = c("Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp",
             "Slp/D", "D/B-H1")
)

p_tem <- tlib |>
  ggplot(aes(x = library, y = prop, fill = temporal_identity)) +
  geom_col() +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(
    values = c("Hth" = "#C76F6B",
               "Hth/Opa" = "#EFC8B9",
               "Opa/Erm" = "#FEE699",
               "Erm/Ey" = "#79A68C",
               "Ey/Hbn" = "#82C9C5",
               "Hbn/Opa/Slp" = "#B4C6E6",
               "Slp/D" = "#A6A1CD",
               "D/B-H1" = "#A46690")
  ) +
  theme_classic() +
  labs(
    x = "Library",
    y = "% of Library",
    fill = "Temporal Identity"
  )

clib <- metadata[, .N, by = c("library", "class")][
  , lib_size := sum(N), by = "library"
][, prop := N / lib_size]

p_class <- clib |>
  ggplot(aes(x = library, y = prop, fill = class)) +
  geom_col() +
  geom_text(aes(label = N), position = position_stack(vjust = 0.5)) +
  scale_fill_manual(
    values = c(
      "NE" = "#E41A1C",
      "NB" = "#377EB8",
      "GMC" = "#4DAF4A",
      "neuron" = "#984EA3",
      "LPC" = "#FF7F00",
      "glia" = "#FFFF33"
    )
  ) +
  theme_classic() +
  labs(
    x = "Library",
    y = "% of Library",
    fill = "Cell Class"
  )

p_ncrna <- metadata |>
  ggplot(aes(x = library, y = nCount_RNA)) +
  geom_violin(aes(fill = library), scale = "area", bounds = c(0, 1e4)) +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "# of RNA UMIs",
    y = "# of RNA UMIs"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )

p_nfrna <- metadata |>
  ggplot(aes(x = library, y = nFeature_RNA)) +
  geom_violin(aes(fill = library), scale = "area") +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "# of RNA Features",
    y = "# of RNA Features"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )

p_ncatac <- metadata |>
  ggplot(aes(x = library, y = nCount_ATAC)) +
  geom_violin(aes(fill = library), scale = "area", bounds = c(0, 2e4)) +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "# of ATAC Insertion Events",
    y = "# of ATAC Insertion Events"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )

p_nfatac <-metadata |>
  ggplot(aes(x = library, y = nFeature_ATAC)) +
  geom_violin(aes(fill = library), scale = "area") +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "# of ATAC Features",
    y = "# of ATAC Features"
  ) +
  theme_classic()+
  theme(
    axis.title.x = element_blank()
  ) +
  theme(
    axis.title.x = element_blank()
  )

p_nuc <- metadata |>
  ggplot(aes(x = library, y = nucleosome_signal)) +
  geom_violin(aes(fill = library), scale = "area") +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "Nucleosome Signal",
    y = "Nucleosome Signal"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )

p_tss <- metadata |>
  ggplot(aes(x = library, y = TSS.enrichment)) +
  geom_violin(aes(fill = library), scale = "area") +
  scale_fill_brewer(palette = "Pastel1") +
  guides(fill = "none") +
  labs(
    title = "TSS Enrichment",
    y = "TSS Enrichment"
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  )


if (!dir.exists("sup_fig/s2")) dir.create("sup_fig/s2", recursive = TRUE)

out1 <- (p_lib / p_spa) | (p_tem / p_class)
ggsave(
  filename = file.path("sup_fig", "s2", "lib_comp.pdf"),
  out1, width = 12.5, height = 10
)

out2 <- (p_ncrna / p_nfrna) | (p_ncatac / p_nfatac) | (p_nuc / p_tss)
ggsave(
  filename = file.path("sup_fig", "s2", "qc_metric.pdf"),
  out2, width = 12, height = 6
)
