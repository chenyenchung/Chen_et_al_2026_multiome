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
  dir.exists(args$proot),
  !is.null(args$depth),
  file.exists(args$depth)
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
  "16" = "LPC",
  "11" = "NE",
  "14" = "glia",
  "41" = "glia",
  "3" = "NB",
  "4" = "NB",
  "0" = "GMC",
  "1" = "GMC",
  "9" = "GMC",
  "6" = "GMC",
  "25" = "neuron",
  "34" = "neuron",
  "10" = "neuron",
  "32" = "neuron",
  "35" = "neuron",
  "40" = "neuron",
  "31" = "neuron",
  "24" = "neuron",
  "39" = "neuron",
  "30" = "neuron",
  "38" = "neuron",
  "42" = "neuron",
  "26" = "neuron",
  "8" = "neuron",
  "36" = "neuron",
  "12" = "neuron",
  "2" = "neuron",
  "22" = "neuron",
  "21" = "neuron",
  "19" = "neuron",
  "13" = "neuron",
  "37" = "neuron",
  "33" = "neuron",
  "23" = "neuron",
  "20" = "neuron",
  "29" = "neuron",
  "5" = "neuron",
  "18" = "neuron",
  "28" = "neuron",
  "7" = "neuron",
  "17" = "neuron",
  "15" = "neuron",
  "27" = "neuron"
)

obj@meta.data$class <- factor(
  class_lut[as.character(Idents(obj))],
  levels = c("NE", "NB", "GMC", "neuron", "LPC", "glia")
)

p_class <- DimRasPlot(obj, group.by = "class") +
  scale_color_manual(
    values = c(
      "NE" = "#1B9E77",
      "NB" = "#D95F02",
      "GMC" = "#7570B3",
      "neuron" = "#66A61E",
      "LPC" = "#E7298A",
      "glia" = "#A6761D"
    )
  ) +
  labs(color = "Cell Class") +
  theme(plot.title = element_blank())

# Notch status: Using Ap as a proxy.
n_lut <- c(
  "16" = NA,
  "11" = NA,
  "14" = NA,
  "41" = NA,
  "3" = NA,
  "4" = NA,
  "0" = NA,
  "1" = NA,
  "9" = NA,
  "6" = NA,
  "25" = "Off",
  "34" = "Off",
  "10" = "On",
  "32" = "Off",
  "35" = "Off",
  "40" = "On",
  "31" = "Off",
  "24" = "On",
  "39" = "On",
  "30" = "On",
  "38" = "Off",
  "42" = "Off",
  "26" = "Off",
  "8" = "On",
  "36" = "On",
  "12" = "On",
  "2" = "On",
  "22" = "Off",
  "21" = "Off",
  "19" = "On",
  "13" = "Off",
  "37" = "Off",
  "33" = "Off",
  "23" = "Off",
  "20" = "On",
  "29" = "On",
  "5" = "On",
  "18" = "On",
  "28" = "On",
  "7" = "Off",
  "17" = "Off",
  "15" = "On",
  "27" = "Off"
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

p_sori <- DimRasPlot(obj, group.by = "spatial_origin") +
  scale_color_manual(
    values = c(
      "pxb" = "#1B9E77",
      "optix" = "#7570B3",
      "dpp" = "#D95F02"
    ),
    na.value = "grey80",
    breaks = c("pxb", "optix", "dpp")
  ) +
  labs(color = "Spatial Origin") +
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
  "16" = NA, # Tll/dac -- LPC
  "11" = NA, # hth/tll/dac -- lateral OPC NE
  "14" = NA, # lamina glia
  "41" = NA, # lamina glia
  "3" = "Erm/Ey", # Opa, erm, (some Hth) but also Hbn and Ey.
  "4" = "Slp/D", # Slp1/Slp2/D
  "0" = "Ey/Hbn", # Hbn and Ey but also Opa and erm.
  "1" = "D/B-H1", # Slp1/Slp2/D/B-H1
  "9" = "Slp/D", # Slp1/Slp2/D/B-H1
  "6" = "Hbn/Opa/Slp", # Ey/Hbn/Opa
  "25" = "Hth/Opa", # Hth/Opa/Run
  "34" = "D/B-H1", # Ets65A, fd95A, oc (T1?)
  "10" = "Slp/D", # Toy, Sox102F, D (Might also have Hbn/Opa/Slp; Tm5Y, cluster 76)
  "32" = "Hth", # Svp, TfAP-2, Hth (Cluster 180?)
  "35" = "Hbn/Opa/Slp", # kn, toy, tj (TmY14, MeSps)
  "40" = "D/B-H1", # Ets65A, TfAP-2 (Cluster 24 / new 50, 75, 76, 77, 256)
  "31" = "Hbn/Opa/Slp", # vvl, tj, fkh, fd59A (Dm8/Dm11/DRA-Dm)
  "24" = "Erm/Ey", # Vvl, Dac (Tm9)
  "39" = "Hth", # Dimm, Bsh, TfAP-2 -- this is probably a mix of Hth - Erm/Ey (TEs)
  "30" = "Hth/Opa", # Erm / TfAP-2 (Tm1/2/4/6)
  "38" = "Hbn/Opa/Slp", # Vvl, tj, fkh, opa, hbn (LPi4-3?)
  "42" = "Hbn/Opa/Slp", # fkh, oc, hbn, fd59A (Dm9, Lai)
  "26" = "Ey/Hbn", # Ey/Tup/Lim3 (Dm10?)
  "8" = "Hbn/Opa/Slp", # Toy, Sox102F, might also have Slp/D? (TmY4, Tm5ab, Tm25...)
  "36" = "Erm/Ey", # Erm, vvl (Tm27(Y), TmY8, Mi10)
  "12" = "Erm/Ey", # Erm, vvl (Tm27(Y), TmY8, Mi10)
  "2" = "Hth/Opa", # Erm / TfAP-2 (Tm1/2/4/6)
  "22" = "Opa/Erm", # Tup, D (without Ey), Lim3 (Mi9)
  "21" = "Hth/Opa", # Run, (also Ey, Tup, D) (Mi4, Pm4)
  "19" = "Hbn/Opa/Slp", # Toy, Sox102F, might also have Slp/D? (TmY4, Tm5ab, Tm25...)
  "13" = "Erm/Ey", # Ey, Kn, Lim3 (Tm29?)
  "37" = "Slp/D", # Hbn, Ets65A, tup (Dm3?)
  "33" = "D/B-H1", # Ets65A, fd95A, Dll (Dm2 and Mi15?)
  "23" = "D/B-H1", # Ets65A, fd95A, oc (T1?)
  "20" = "Hth", # Hth, Bsh (Mi1)
  "29" = "Hth/Opa", # Hth/Run/Bsh
  "5" = "Hbn/Opa/Slp", # vvl, tj, fkh, fd59A (Dm8/Dm11/DRA-Dm)
  "18" = "Opa/Erm", # Tup, D (without Ey), Lim3 (Mi9) (Could have Mi4/Pm4 as well)
  "28" = "D/B-H1", # Ets65A/Dll/Svp
  "7" = "Slp/D", # Dll, fd59A, Ets65A (Cluster 138?)
  "17" = "Hbn/Opa/Slp", # vvl, tj, fkh, fd59A (Dm8/Dm11/DRA-Dm)
  "15" = "Hbn/Opa/Slp", # toy/Sox102F
  "27" = "D/B-H1" # fd59A, oc, B-H1
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
ggsave_sep(p_nstat, file.path("fig", "notch_status.pdf"), width = 5, height = 4)
ggsave_sep(p_tid, file.path("fig", "temporal_identity.pdf"), width = 5, height = 4)
ggsave_sep(p_sori, file.path("fig", "spatial_origin.pdf"), width = 5, height = 4)
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

idepth <- read.csv(args$depth, header = FALSE, row.names = 1)
colnames(idepth) <- "nCount_insertion"
obj <- AddMetaData(obj, metadata = idepth)
obj@meta.data$X <- NULL

# Compute mixing statistics before saving
obj[["seurat_clusters"]] <- Idents(obj)
Idents(obj) <- "class"

get_global_simpsons <- function(obj, group.by) {
  label_freq <- table(obj[[]][[group.by]])
  label_proportion <- label_freq / length(obj[[]][[group.by]])
  return(sum(label_proportion ^ 2))
}

get_local_simpsons <- function(obj, group.by, nn.name, gsi) {
  lsi <- vapply(seq_len(ncol(obj)), function(x) {
    idx <- obj@neighbors[[nn.name]]@nn.idx[x, ]
    local_freq <- table(obj[[]][[group.by]][idx])
    local_prop <- local_freq / length(idx)
    local_si <- sum(local_prop ^ 2)
    return((local_si - gsi) / (1 - gsi))
  }, numeric(1))
}

spatial_gsi <- get_global_simpsons(obj, "spatial_origin")
spatial_clsi <- get_local_simpsons(
  obj, "spatial_origin", "weighted.nn", spatial_gsi
)

# Remove NE and LPC (no temporal annotation)
temporal_gsi <- get_global_simpsons(
  subset(obj, idents = c("NB", "GMC", "neuron")), "temporal_identity"
)
temporal_clsi <- get_local_simpsons(
  obj, "temporal_identity", "weighted.nn", temporal_gsi
)

temporal_clsi[!Idents(obj) %in% c("NB", "GMC", "neuron")] <- NA

# Only neurons has Notch annotation
notch_gsi <- get_global_simpsons(
  subset(obj, idents = "neuron"), "notch_status"
)

notch_clsi <- get_local_simpsons(
  obj, "notch_status", "weighted.nn", notch_gsi
)
notch_clsi[!Idents(obj) == "neuron"] <- NA

obj$spatial_clsi <- spatial_clsi
obj$temporal_clsi <- temporal_clsi
obj$notch_clsi <- notch_clsi

Idents(obj) <- "seurat_clusters"

saveRDS(obj, file.path("obj", "annotated.rds"))
write.csv(obj@meta.data, "metadata.csv")
