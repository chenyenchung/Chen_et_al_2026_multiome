#!/usr/bin/env Rscript
library(R.utils)
options(future.globals.maxSize = 1024 ^ 3)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$obj <-"int/filtered_obj.rds"
  args$p <- "int/hvg_np/varp.rds"
  args$g <- "int/hvg_np/varg.rds"
  args$proot <- "/scratch/ycc520/thesis"
  args$resolution <- "0.5"
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
  !is.na(as.numeric(args$resolution))
)

renv::load(args$proot)

library(Seurat)
library(Signac)
library(harmony)

obj <- readRDS(args$obj)
varp <- readRDS(args$p)
varg <- readRDS(args$g)

p_tbl <- table(unlist(varp))
g_tbl <- table(unlist(varg))

g <- names(g_tbl[g_tbl > 2])
p <- names(p_tbl[p_tbl > 2])

DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj,
                   vars.to.regress = "percent.mt",
                   residual.features = g,
                   return.only.var.genes = FALSE)
obj <- RunPCA(obj, npcs = 400, assay = "SCT")
obj <- RunHarmony(
  obj, group.by.vars = "library", reduction.use = "pca", reduction.save = "hpca"
)

DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- RunSVD(obj, features = p, n = 400, assay = "ATAC")
obj <- RunHarmony(
  obj, group.by.vars = "library", reduction.use = "lsi",
  reduction.save = "hlsi", project.dim = FALSE
)

obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("hpca", "hlsi"), dims.list = list(1:250, 1:250)
)
obj <- FindClusters(
  obj, resolution = as.numeric(args$resolution), graph.name = "wsnn"
)
obj <- RunUMAP(obj, nn.name = "weighted.nn")

saveRDS(obj, "filtered_wnn.rds")
