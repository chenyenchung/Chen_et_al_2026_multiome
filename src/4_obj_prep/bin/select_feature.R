#!/usr/bin/env Rscript
FindVariablePeaks <- function(obj, n, assay = "ATAC") {
  tfidf_mat <- LayerData(obj, layer = "data", assay = assay)
  peak_var <- Matrix::rowMeans(tfidf_mat^2) - Matrix::rowMeans(tfidf_mat) ^ 2
  out <- sort(peak_var, decreasing = TRUE) |> head(n = n)
  return(names(out))
}

library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "/scratch/ycc520/thesis/int/base_obj.rds"
  args$iter <- "5"
  args$resolution <- "0.2"
}

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$iter),
  !is.na(as.integer(args$iter)),
  !is.null(args$resolution),
  !is.na(as.numeric(args$resolution))
)

renv::load(args$proot)

library(Seurat)
library(Signac)
library(Matrix)

obj <- readRDS(args$obj)
obj <- subset(obj, permissive_filter)

DefaultAssay(obj) <- "RNA"
obj <- PercentageFeatureSet(obj, pattern = "^mt:", col.name = "percent.mt")
obj <- SCTransform(obj, vars.to.regress = "percent.mt")
obj <- RunPCA(obj, npcs = 100)
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj, n = 101)

obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:100, 2:101))
obj <- FindClusters(obj, resolution = as.numeric(args$resolution), graph.name = "wsnn")
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "umap.initial")

ovarg <- VariableFeatures(obj, assay = "SCT")
ovarp <- FindVariablePeaks(obj, n = 10000)

var_features <- list(g = list(), p = list())
max_iter <- as.integer(args$iter)

for (i in seq_len(max_iter)) {
  obj_list <- SplitObject(obj, split.by = paste0("wsnn_res.", args$resolution))
  varl <- lapply(
    obj_list, function(x) {
      x <- SCTransform(x, vars.to.regress = "percent.mt", return.only.var.genes = FALSE)
      p <- FindVariablePeaks(x, n = 3000)
      Idents(x) <- "spatial_origin"
      
      scontrasts <- list(c("pxb", "optix"), c("optix", "dpp"), c("pxb", "dpp"))
      sg <- lapply(scontrasts, function(contrast) {
        rmarkers <- FindMarkers(
          x, ident.1 = contrast[1], ident.2 = contrast[2], assay = "SCT"
        )
        sig_m <- subset(rmarkers, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
        return(row.names(sig_m))
      }
      )
      
      sp <- lapply(scontrasts, function(contrast) {
        rmarkers <- FindMarkers(
          x, ident.1 = contrast[1], ident.2 = contrast[2], assay = "ATAC"
        )
        sig_m <- subset(rmarkers, p_val_adj < 0.05 & abs(avg_log2FC) > 1)
        return(row.names(sig_m))
      }
      )
      
      sgout <- unique(c(sg, VariableFeatures(x)))
      spout <- unique(c(sp, p))
      return(list(g = sgout, p = spout))
    }
  )
  rm(obj_list)
  gc()
  
  var_features$g[[i]] <- lapply(varl, function(x) return(x$g)) |>
    unlist() |> unique()
  var_features$p[[i]] <- lapply(varl, function(x) return(x$p)) |>
    unlist() |> unique()
  
  vg_tbl <- table(unlist(var_features$g))
  vp_tbl <- table(unlist(var_features$p))
  
  g <-  names(vg_tbl)
  p <-  names(vp_tbl)
  

  obj <- SCTransform(obj, vars.to.regress = "percent.mt",
                     residual.features = g)
  obj <- RunPCA(obj, npcs = 100, assay = "SCT")
  obj <- RunSVD(obj, features = p, n = 101, assay = "ATAC")
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:100, 2:101))
  obj <- FindClusters(obj, resolution = as.numeric(args$resolution), graph.name = "wsnn")
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = paste0("umap_", i))
}

saveRDS(var_features$g, "varg.rds")
saveRDS(var_features$p, "varp.rds")

