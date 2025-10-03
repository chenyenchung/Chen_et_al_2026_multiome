#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$obj <-"int/base_obj.rds"
  args$p <- "work/cd/be383e1a1285346200bb6c7673d158/varp.rds"
  args$g <- "work/cd/be383e1a1285346200bb6c7673d158/varg.rds"
  args$proot <- "/scratch/ycc520/thesis"
  args$resolution <- "1"
  args$nnout <- "int/nn_pred"
  args$anno <- "static/ozel_2020_lut.csv"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$p),
  file.exists(args$p),
  !is.null(args$g),
  file.exists(args$g),
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$resolution),
  !is.na(as.numeric(args$resolution)),
  !is.null(args$nnout),
  dir.exists(args$nnout),
  !is.null(args$anno),
  file.exists(args$anno)
)

renv::load(args$proot)

library(Seurat)
library(Signac)

obj <- readRDS(args$obj)
varp <- readRDS(args$p)
varg <- readRDS(args$g)

p_tbl <- table(unlist(varp))
g_tbl <- table(unlist(varg))

g <- names(g_tbl[g_tbl > 2])
p <- names(p_tbl[p_tbl > 2])

obj <- PercentageFeatureSet(obj, pattern = "^mt:", col.name = "percent.mt")

# Get label-annotation conversion
anno_lut <- read.csv(args$anno)
anno_lut$cluster <- as.character(anno_lut$cluster)
lut <- anno_lut$anno
names(lut) <- anno_lut$cluster

obj <- AddMetaData(
  obj,
  metadata = readLines(file.path(args$nnout, "TargetPreds.txt")),
  col.name = "nn_pred_id"
)
obj <- AddMetaData(
  obj,
  metadata = unname(lut[readLines(file.path(args$nnout, "TargetPreds.txt"))]),
  col.name = "nn_pred_label"
)
obj <- AddMetaData(
  obj,
  metadata = as.numeric(readLines(file.path(args$nnout, "TargetPreds-Confidence.txt"))),
  col.name = "nn_conf"
)

obj <- subset(obj, permissive_filter)

DefaultAssay(obj) <- "RNA"

obj <- SCTransform(obj,
                   vars.to.regress = "percent.mt",
                   residual.features = g,
                   return.only.var.genes = FALSE)
obj <- RunPCA(obj, npcs = 100, assay = "SCT")
obj <- RunSVD(obj, features = p, n = 101, assay = "ATAC")
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:100, 2:101))
obj <- FindClusters(obj, resolution = as.numeric(args$resolution), graph.name = "wsnn")
obj <- RunUMAP(obj, nn.name = "weighted.nn")

saveRDS(obj, "base_ann_obj.rds")
