#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <-"int/raw_cluster.rds"
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
library(ggrastr)
library(data.table)
library(patchwork)

DimRasPlot <- function(
    object, reduction = "umap", group.by = NULL, pt.size = 0.2) {
  if (is.null(group.by)) {
    group.by = "seurat_clusters"
  }
  
  plot_df <- Embeddings(object@reductions[[reduction]]) |>
    as.data.frame()
  colnames(plot_df) <- paste0("dim_", seq_len(ncol(plot_df)))
  
  plot_df[[group.by]] <- object[[group.by]][, 1]
  p <- plot_df |>
    ggplot(aes(x = dim_1, y = dim_2, color = .data[[group.by]])) +
    geom_point_rast(size = pt.size) +
    theme_void()
  return(p)
}
ggsave_sep <- function(
    plot, filename,
    width = NULL, height = NULL, limitsize = FALSE
    ) {
  legendsp <- cowplot::get_plot_component(plot, "guide-box", return_all = TRUE)
  plot <- plot + theme(legend.position="none")
  
  if (!inherits(legendsp, "gtable")) {
    ## Drop empty elements
    to_keep <- sapply(legendsp, function(x) !"zeroGrob" %in% class(x))
    legendsp <- legendsp[to_keep]
    if (length(legendsp) == 1) {
      legendsp <- legendsp[[1]]
    } else {
      stop("Legend extraction error: There is more than 1 item.")
    }
  }
  basename <- strsplit(filename, "\\.")[[1]][1]
  ext <- strsplit(filename, "\\.")[[1]][2]
  ggsave(
    filename = filename, plot = plot,
    width = width , height = height, limitsize = limitsize
  )
  ggsave(
    filename = paste0(basename, "_guide.", ext),
    legendsp
  )
}

obj <- readRDS(args$obj)

# Class annotation
# NE: Grh and shg
# NB: Mira
# GMC: Ase
# Neuron: CadN
# Glia: Repo and AnxB9

class_plot <- AggregateExpression(
  obj, assay = "RNA", features = c(
    "grh", "shg", "mira", "ase", "N", "CadN", "repo", "AnxB9", "gcm"
  ),
  return.seurat = TRUE
)
class_plot <- NormalizeData(class_plot)
class_plot <- ScaleData(class_plot)
class_dt <- data.frame(LayerData(class_plot, layer = "scale.data"))
class_dt$gene <- row.names(class_dt)
class_dt <- as.data.table(class_dt)
class_dt <- melt(
  class_dt, value.name = "z_score", id.vars = "gene", variable.name = "cluster"
)

cluster_order <- t(LayerData(class_plot, layer = "scale.data")) |>
  dist() |> hclust()
gene_order <- LayerData(class_plot, layer = "scale.data") |>
  dist() |> hclust()

class_mk_hm <- class_dt[
  , cluster := factor(
    cluster,
    levels = cluster_order$labels[cluster_order$order],
    labels = sub("^g", "", cluster_order$labels[cluster_order$order])
  )
][
  , gene := factor(
    gene,
    levels = gene_order$labels[gene_order$order]
  )
] |>
  ggplot(aes(x = cluster, y = gene, fill = z_score)) +
  geom_tile() +
  labs(x = "Cluster", fill = "Scaled Exp.\n(Z-score)") +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  scale_fill_gradient2(low = "#00441B", mid = "#F7F7F7", high = "#40004B")

class_lut <- c(
  "11" = "glia",
  "8" = "NE",
  "0" = "GMC",
  "4" = "NB",
  "15" = "glia",
  "10" = "neuron",
  "12" = "neuron",
  "7" = "neuron",
  "1" = "neuron",
  "9" = "neuron",
  "3" = "neuron",
  "6" = "neuron",
  "2" = "neuron",
  "5" = "neuron",
  "13" = "neuron",
  "14" = "neuron"
)

obj@meta.data$class <- factor(
  class_lut[as.character(Idents(obj))],
  levels = c("NE", "NB", "GMC", "neuron", "glia")
)

p_class <- DimRasPlot(obj, group.by = "class") +
  scale_color_manual(
    values = c(
      "NE" = "#1B9E77",
      "NB" = "#D95F02",
      "GMC" = "#7570B3",
      "neuron" = "#66A61E",
      "glia" = "#E7298A"
    )
  ) +
  labs(color = "Cell Class") +
  theme(plot.title = element_blank())

bclass_lut <- c(
  "11" = "glia",
  "8" = "progenitor",
  "0" = "progenitor",
  "4" = "progenitor",
  "15" = "glia",
  "10" = "neuron",
  "12" = "neuron",
  "7" = "neuron",
  "1" = "neuron",
  "9" = "neuron",
  "3" = "neuron",
  "6" = "neuron",
  "2" = "neuron",
  "5" = "neuron",
  "13" = "neuron",
  "14" = "neuron"
)
obj@meta.data$bclass <- factor(
  bclass_lut[as.character(Idents(obj))],
  levels = c("progenitor", "neuron", "glia")
)

p_bclass <- DimRasPlot(obj, group.by = "bclass") +
  scale_color_manual(
    values = c(
      "progenitor" = "#E6AB02",
      "neuron" = "#66A61E",
      "glia" = "#E7298A"
    )
  ) +
  labs(color = "Broad Class") +
  theme(plot.title = element_blank())

n_lut <- c(
  "11" = NA,
  "8" = NA,
  "0" = NA,
  "4" = NA,
  "15" = NA,
  "10" = "Notch_OFF",
  "12" = "Notch_OFF",
  "7" = "Notch_OFF",
  "1" = "Notch_ON",
  "9" = "Notch_OFF",
  "3" = "Notch_ON",
  "6" = "Notch_OFF",
  "2" = "Notch_ON",
  "5" = "Notch_ON",
  "13" = "Notch_ON",
  "14" = "Notch_OFF"
)

obj@meta.data$notch_status <- factor(
  n_lut[as.character(Idents(obj))],
  levels = c("Notch_OFF", "Notch_ON")
)

p_nstat <- DimRasPlot(obj, group.by = "notch_status") +
  scale_color_manual(
    values = c(
      "Notch_OFF" = "#435274",
      "Notch_ON" = "#ba3c3c"
    ),
    na.value = "grey80",
    breaks = c("Notch_OFF", "Notch_ON"),
    labels = c("Off", "On")
  ) +
  labs(color = "Notch Status") +
  theme(plot.title = element_blank())

on_obj <- subset(obj, idents = c(
  "11", "8", "0", "4", "15", "1", "3", "2", "5", "13"
))

off_obj <- subset(obj, idents = c(
  "11", "8", "0", "4", "15", "10", "12", "7", "9", "6", "14"
))


# Output
if (!dir.exists("sup")) {
  dir.create("sup")
}
if (!dir.exists("fig")) {
  dir.create("fig")
}
if (!dir.exists("obj")) {
  dir.create("obj")
}
ggsave_sep(p_class, file.path("fig", "class.pdf"), width = 5, height = 4)
ggsave_sep(p_bclass, file.path("sup", "broad_class.pdf"), width = 5, height = 4)
ggsave_sep(p_nstat, file.path("fig", "notch_status.pdf"), width = 5, height = 4)
ggsave_sep(
  class_mk_hm, file.path("sup", "class_marker_heatmap.pdf"),
  width = 6, height = 4
)
saveRDS(obj, file.path("obj", "obj_cl_n.rds"))
saveRDS(on_obj, file.path("obj", "non_obj.rds"))
saveRDS(off_obj, file.path("obj", "noff_obj.rds"))
