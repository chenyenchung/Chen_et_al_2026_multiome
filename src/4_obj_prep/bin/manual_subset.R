#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
if (interactive()) {
  args$obj <-"int/objects/base_ann_obj.rds"
  args$proot <- "/scratch/ycc520/thesis"
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
library(AnnotationHub)
library(patchwork)


obj_path <- args$obj
obj <- readRDS(obj_path)

# DimPlot(obj, label = TRUE) + NoLegend() |
#   VlnPlot(obj, "nn_conf", pt.size = 0) + NoLegend()

anno_update <- as.POSIXct("2026-01-12 00:55:00 EST")
obj_creation <- file.info(obj_path)$ctime

# If the object is newer than annotation, we should manually review
stopifnot(
  obj_creation < anno_update
)

# table(obj$wsnn_res.1, obj$nn_pred_label)["27", ] |> sort()
# a <- FindMarkers(obj, "29", only.pos = TRUE, assay = "SCT")
# a <- subset(a, p_val_adj < 0.001 & avg_log2FC > 1.5)

exclude <- c(
  # Low confidence overall
  "21", # Closest call: 85
  "11", # Closest call: 85
  "36", # Closest call: 102
  "18", # Closest call: 85
  "6", # Closest call: 102
  "30", # Closest call: 85/102
  "37", # Unclear cells expressing Ilp3, wrapping glia?
  "13", # Ase & Run+ & Optix+ -- likely not OPC-origin
  "28", # Likely LC10b
  "12", # No positive marker at adj.p < 0.001 & log2FC > 1.5, suspecting LQ.
  "26", # IPC progenitors (GMC1/GMC2 like)?
  "35", # Perineural glia?
  "20", # Unknown Obp-rich population -- likely not OPC-origin.
  "33", # Likely C2/3
  "34", # Likely T2/3
  "27" # Likely T4/5
)

obj <- subset(obj, idents = setdiff(levels(Idents(obj)), exclude))

DefaultAssay(obj) <- "RNA"
obj <- SCTransform(
  obj, vars.to.regress = "percent.mt", return.only.var.genes = FALSE
)
obj <- RunPCA(obj, npcs = 100)
DefaultAssay(obj) <- "ATAC"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj, n = 101)

obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("pca", "lsi"), dims.list = list(1:100, 2:101)
)
obj <- RunUMAP(obj, nn.name = "weighted.nn")

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

cbc_list <- split(row.names(obj@meta.data), obj$library)

f_list <- lapply(names(ins_list), function(lib) {
  f <- CreateFragmentObject(
    file.path("./", ins_list[[lib]]),
    cells = cbc_list[[lib]]
  )
  return(f)
})

DefaultAssay(obj) <- "ATAC"
Fragments(obj) <- f_list
obj <- NucleosomeSignal(object = obj)
obj <- TSSEnrichment(object = obj)
obj <- subset(obj, TSS.enrichment > 2 & obj$nucleosome_signal < 4)

write.csv(obj@meta.data, "metadata.csv")

saveRDS(obj, "intermediate.rds")
