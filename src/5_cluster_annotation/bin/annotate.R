#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <-"int/objects/raw_cluster.rds"
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
library(AnnotationHub)

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
    guides(color = guide_legend(override.aes = list(size = 6))) +
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
    "dac", "eya" ,"grh", "shg", "mira", "ase", "pros", "insc", "N", "CadN",
    "repo", "gcm"
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
  "20" = "LPC",
  "15" = "NE",
  "2" = "NB",
  "4" = "NB",
  "0" = "GMC",
  "1" = "GMC",
  "9" = "GMC",
  "5" = "GMC",
  "17" = "GMC",
  "30" = "GMC",
  "10" = "neuron",
  "3" = "neuron",
  "34" = "neuron",
  "19" = "neuron",
  "31" = "neuron",
  "27" = "neuron",
  "25" = "neuron",
  "12" = "neuron",
  "8" = "neuron",
  "7" = "neuron",
  "14" = "neuron",
  "26" = "neuron",
  "32" = "neuron",
  "23" = "neuron",
  "6" = "neuron",
  "33" = "neuron",
  "18" = "neuron",
  "21" = "neuron",
  "11" = "neuron",
  "24" = "neuron",
  "16" = "neuron",
  "22" = "neuron",
  "28" = "neuron",
  "29" = "neuron",
  "35" = "neuron",
  "13" = "neuron"
)

obj@meta.data$class <- factor(
  class_lut[as.character(Idents(obj))],
  levels = c("NE", "NB", "GMC", "neuron", "LPC")
)

p_class <- DimRasPlot(obj, group.by = "class") +
  scale_color_manual(
    values = c(
      "NE" = "#1B9E77",
      "NB" = "#D95F02",
      "GMC" = "#7570B3",
      "neuron" = "#66A61E",
      "LPC" = "#E7298A"
    )
  ) +
  labs(color = "Cell Class") +
  theme(plot.title = element_blank())

# Broad classes: Pooling NE, NB, and GMC
bclass_lut <- c(
  "20" = "LPC",
  "15" = "progenitor",
  "2" = "progenitor",
  "4" = "progenitor",
  "0" = "progenitor",
  "1" = "progenitor",
  "9" = "progenitor",
  "5" = "progenitor",
  "17" = "progenitor",
  "30" = "progenitor",
  "3" = "neuron",
  "10" = "neuron",
  "34" = "neuron",
  "19" = "neuron",
  "31" = "neuron",
  "27" = "neuron",
  "25" = "neuron",
  "12" = "neuron",
  "8" = "neuron",
  "7" = "neuron",
  "14" = "neuron",
  "26" = "neuron",
  "32" = "neuron",
  "23" = "neuron",
  "6" = "neuron",
  "33" = "neuron",
  "18" = "neuron",
  "21" = "neuron",
  "11" = "neuron",
  "24" = "neuron",
  "16" = "neuron",
  "22" = "neuron",
  "28" = "neuron",
  "29" = "neuron",
  "35" = "neuron",
  "13" = "neuron"
)
obj@meta.data$bclass <- factor(
  bclass_lut[as.character(Idents(obj))],
  levels = c("progenitor", "neuron", "LPC")
)

p_bclass <- DimRasPlot(obj, group.by = "bclass") +
  scale_color_manual(
    values = c(
      "progenitor" = "#E6AB02",
      "neuron" = "#66A61E",
      "LPC" = "#E7298A"
    )
  ) +
  labs(color = "Broad Class") +
  theme(plot.title = element_blank())

# Notch status: Using Ap as a proxy, but cluster 24 is an exception.
# Cite Coyne et al. 2025.
n_lut <- c(
  "20" = NA,
  "15" = NA,
  "2" = NA,
  "4" = NA,
  "0" = NA,
  "1" = NA,
  "9" = NA,
  "5" = NA,
  "17" = NA,
  "30" = NA,
  "10" = "On",
  "3" = "On",
  "34" = "Off",
  "19" = "Off",
  "31" = "Off",
  "27" = "Off",
  "25" = "On",
  "12" = "On",
  "8" = "On",
  "7" = "On",
  "14" = "On",
  "26" = "Off",
  "32" = "Off",
  "23" = "On",
  "6" = "On",
  "33" = "On",
  "18" = "Off",
  "21" = "Off",
  "11" = "On",
  "24" = "On",
  "16" = "Off",
  "22" = "Off",
  "28" = "Off",
  "29" = "Off",
  "35" = "Off",
  "13" = "Off"
)

obj@meta.data$notch_status <- factor(
  n_lut[as.character(Idents(obj))],
  levels = c("Off", "On")
)

p_nstat <- DimRasPlot(obj, group.by = "notch_status") +
  scale_color_manual(
    values = c(
      "Off" = "#435274",
      "On" = "#ba3c3c"
    ),
    na.value = "grey80",
    breaks = c("Off", "On")
  ) +
  labs(color = "Notch Status") +
  theme(plot.title = element_blank())

# Temporal identity: By concentric genes
conc_plot <- AggregateExpression(
  obj, assay = "RNA", features = c(
    "hth", "svp", "bsh", "TfAP-2", "run", "Lim3", "tup", "vvl", "kn", "toy",
    "Sox102F", "fd59A", "tj", "fkh", "oc", "hbn", "Ets65A", "Dll", "erm",
    "dimm", "opa", "ey", "slp1", "slp2", "D", "B-H1", "tll", "dac"
  ),
  return.seurat = TRUE
)
conc_plot <- NormalizeData(conc_plot)
conc_plot <- ScaleData(conc_plot)
conc_dt <- data.frame(LayerData(conc_plot, layer = "scale.data"))
conc_dt$gene <- row.names(conc_dt)
conc_dt <- as.data.table(conc_dt)
conc_dt <- melt(
  conc_dt, value.name = "z_score", id.vars = "gene", variable.name = "cluster"
)

cluster_order <- t(LayerData(conc_plot, layer = "scale.data")) |>
  dist() |> hclust()
gene_order <- LayerData(conc_plot, layer = "scale.data") |>
  dist() |> hclust()

conc_mk_hm <- conc_dt[
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
][
  , z_score := pmin(z_score, 1.5)
][
  , z_score := pmax(z_score, -1.5)
] |>
  ggplot(aes(x = cluster, y = gene, fill = z_score)) +
  geom_tile() +
  labs(x = "Cluster", fill = "Scaled Exp.\n(Z-score)") +
  theme_classic() +
  theme(axis.title.y = element_blank()) +
  scale_fill_gradient2(
    low = "#00441B", mid = "#F7F7F7", high = "#40004B"
  )

t_lut <- c(
  "20" = NA, # Tll -- LPC
  "15" = NA, # Hth, Tll -- Lateral OPC NE
  "2" = "D/B-H1", # D, B-H1, Tll
  "4" = "Opa/Erm", # Opa, erm, (some Hth) but also Hbn and Ey.
  "0" = "Slp/D", # Ets65A, svp, fd59A
  "1" = "Erm/Ey", # Ey, erm (Also Opa), but also slp
  "9" = "Hbn/Opa/Slp", # Slp, opa
  "5" = "Ey/Hbn", # Ey, Hbn, but also have Erm
  "17" = "Hbn/Opa/Slp", # Vvl, tj, fkh, opa, hbn
  "30" = "D/B-H1", # Svp, B-H1
  "10" = "Hth", # Hth
  "3" = "Hth/Opa", # Svp, Erm (Also Hbn/Opa/Slp -- likely Tm1/4/2/6)
  "34" = "Opa/Erm", # Tup, D (without Ey), Lim3
  "19" = "Slp/D", # Ets65A, fd95A, Dll, Svp (Cluster 138/ (new 212)?)
  "31" = "Slp/D", # Hbn, Ets65A, tup (Dm3?)
  "27" = "D/B-H1", # Ets65A, fd95A, Dll (Dm2 and Mi15?)
  "25" = "D/B-H1", # Ets65A, TfAP-2 (Cluster 24 / new 50, 75, 76, 77, 256)
  "12" = "Slp/D", # Svp, Ets65A
  "8" = "Hbn/Opa/Slp", # Toy, Sox102F, might also have Slp/D?
  "7" = "Hbn/Opa/Slp", # Toy, Sox102F, might also have Slp/D? (TmY4, Tm5ab, Tm25...)
  "14" = "Slp/D", # Toy, Sox102F, D (Might also have Hbn/Opa/Slp; Tm5Y, cluster 76)
  "26" = "Hbn/Opa/Slp", # vvl, tj, fkh, fd59A (Dm8/Dm11/DRA-Dm)
  "32" = "Hbn/Opa/Slp", # fkh, hbn, tj (Lpi4-3)
  "23" = "Hth", # Hth, Bsh (Mi1)
  "6" = "Hth/Opa", # Erm / TfAP-2 (Tm1/2/4/6)
  "33" = "Hth", # Dimm, Bsh, TfAP-2 -- this is probably a mix of Hth - Erm/Ey (TEs)
  "18" = "Erm/Ey", # Ey, Kn, Lim3 (Tm29?)
  "21" = "Opa/Erm", # Lim3, Tup, D (No Ey & Kn) (Mi9)
  "11" = "Erm/Ey", # Erm, vvl (Tm27(Y), TmY8, Mi10)
  "24" = "Erm/Ey", # Vvl, Dac (Tm9)
  "16" = "Hth/Opa", # Run, (also Ey, Tup, D) (Mi4, Pm4)
  "22" = "Ey/Hbn", # tup, ey, toy (Cluster 173, 181, Dm10)
  "28" = "Hth", # Svp, TfAP-2, Hth (Cluster 180?)
  "29" = "Hbn/Opa/Slp", # kn, toy, tj (TmY14, MeSps)
  "35" = "Hbn/Opa/Slp", # fkh, oc, hbn (Dm9, Lai)
  "13" = "D/B-H1" # Ets65A, fd95A, oc (T1?)
)

obj@meta.data$temporal_identity <- factor(
  t_lut[as.character(Idents(obj))],
  levels = c(
    "Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp",
    "Slp/D", "D/B-H1"
  )
)

p_tid <- DimRasPlot(obj, group.by = "temporal_identity") +
  scale_color_manual(
    values = c("Hth" = "#C76F6B",
               "Hth/Opa" = "#EFC8B9",
               "Opa/Erm" = "#FEE699",
               "Erm/Ey" = "#79A68C",
               "Ey/Hbn" = "#82C9C5",
               "Hbn/Opa/Slp" = "#B4C6E6",
               "Slp/D" = "#A6A1CD",
               "D/B-H1" = "#A46690"),
    na.value = "grey80",
    breaks = c(
      "Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", "Ey/Hbn", "Hbn/Opa/Slp",
      "Slp/D", "D/B-H1"
    )
  ) +
  labs(color = "Temporal Identity") +
  theme(plot.title = element_blank())

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
ggsave_sep(p_tid, file.path("fig", "temporal_identity.pdf"), width = 5, height = 4)
ggsave_sep(
  class_mk_hm, file.path("sup", "class_marker_heatmap.pdf"),
  width = 6, height = 4
)
ggsave_sep(
  conc_mk_hm, file.path("sup", "temporal_marker_heatmap.pdf"),
  width = 6, height = 4
)

ah <- AnnotationHub()
# query(ah, "EnsDb.Dmelanogaster.v88")
# names(): AH53702
# $dataprovider: Ensembl
# $species: Drosophila melanogaster
annotations <- GetGRangesFromEnsDb(ensdb = ah[["AH53702"]])
genome(annotations) <- "dm6"
Annotation(obj) <- annotations

ins_list <- list.files("./", pattern = "\\.tsv\\.gz$")
names(ins_list) <- vapply(ins_list, function(x) {
  return(sub("_all\\.tsv\\.gz$", "", x))
}, FUN.VALUE = character(1))



f_list <- lapply(names(ins_list), function(lib) {
  f <- CreateFragmentObject(
    file.path("./", ins_list[[lib]]),
    cells = row.names(obj@meta.data)[obj$library == lib]
  )
  return(f)
})

DefaultAssay(obj) <- "ATAC"
Fragments(obj) <- f_list

saveRDS(obj, file.path("obj", "annotated.rds"))

on_obj <- subset(obj, idents = c(
  "15", "2", "4", "0", "1", "9", "5", "17", "30",
  "10", "3", "25", "12", "8", "7", "14", "23", "6", "33", "11", "24"
))
saveRDS(on_obj, file.path("obj", "non_obj.rds"))
rm(on_obj)
gc()

off_obj <- subset(obj, idents = c(
  "15", "2", "4", "0", "1", "9", "5", "17", "30",
  "34", "19", "31", "27", "26", "32", "18", "21",
  "16", "22", "28", "29", "35", "13"
))
saveRDS(off_obj, file.path("obj", "noff_obj.rds"))
