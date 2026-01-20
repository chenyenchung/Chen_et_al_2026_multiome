#!/usr/bin/env Rscript
ReadH5mat <- function(libs, bclist, h5Files, ftype = NULL) {
  mlist <- lapply(
    libs, function(sn) {
      # Open H5 file with error handling
      h5fp <- tryCatch(
        H5File$new(filename = h5Files[[sn]], mode = "r"),
        error = function(e) {
          stop(paste("Failed to open H5 file for", sn, ":", conditionMessage(e)))
        }
      )
      
      # Verify expected datasets exist
      required_paths <- c("matrix/data", "matrix/indices", "matrix/indptr",
                          "matrix/shape", "matrix/features/name",
                          "matrix/features/feature_type", "matrix/barcodes")
      for (path in required_paths) {
        if (!h5fp$exists(path)) {
          stop(paste("Missing dataset", path, "in", h5Files[[sn]]))
        }
      }
      
      if (is.null(ftype)) {
        ftype <- h5fp[["matrix/features/feature_type"]][]
        message("Available ftype(s): ", paste(unique(ftype), collapse = " "))
        stop("ftype must be set.")
      }
      
      # From Seurat::Read10X_h5
      counts <- h5fp[["matrix/data"]][]
      indices <- h5fp[["matrix/indices"]][]
      indptr <- h5fp[["matrix/indptr"]][]
      shp <- h5fp[["matrix/shape"]][]
      if (ftype == "Peaks") {
        features <- h5fp[["matrix/features/interval"]][]
      } else {
        features <- h5fp[["matrix/features/name"]][]
      }
      
      ftype <- h5fp[["matrix/features/feature_type"]][]
      barcodes <- h5fp[["matrix/barcodes"]][]
      sparse.mat <- sparseMatrix(
        # Convert row indices: 0-based (HDF5) → 1-based (R)
        i = indices + 1,
        # Column pointers stay 0-based (CSC format requirement)
        p = indptr,
        x = counts,
        dims = shp,
        repr = "C"
      )
      rownames(x = sparse.mat) <- features
      colnames(x = sparse.mat) <- barcodes
      
      
      cbc_pass <- bclist[[sn]]
      cbc_pass <- cbc_pass[cbc_pass %in% colnames(sparse.mat)]
      
      # Validate barcode whitelist is not empty
      if (length(cbc_pass) == 0) {
        stop(paste("No valid barcodes for library", sn))
      }
      
      # Filter features to keep only desired assay (do this before cell filtering)
      # We are rearranging row order later, so converting the matrix to a
      # row-major format will make it faster
      feature_idx <- which(ftype == ftype)
      sparse.mat <- as(sparse.mat, "RsparseMatrix")
      sparse.mat <- sparse.mat[feature_idx, cbc_pass]
      colnames(sparse.mat) <- paste(sn, colnames(sparse.mat), sep = "#")
      
      return(sparse.mat)
    }
  )
  return(mlist)
}

library(R.utils)
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$gtf <- "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf"
  args$data <- "/scratch/ycc520/thesis/data/bam_new"
}

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$data),
  dir.exists(args$data),
  !is.null(args$gtf),
  file.exists(args$gtf)
)

renv::load(args$proot)
# # library(Seurat)
library(Matrix)
# # library(Signac)
library(rtracklayer)
# library(GenomicRanges)
library(GenomeInfoDb)
# library(BSgenome.Dmelanogaster.BDGP.dm6)
# library(data.table)
library(hdf5r)

libs <- c("stf_2", "stf_3", "stf_4", "stf_5")
names(libs) <- libs

## Load permissive barcode lists
bcFiles <- c(
  stf_2 = file.path("filtered_barcode", "stf_2.txt"),
  stf_3 = file.path("filtered_barcode", "stf_3.txt"),
  stf_4 = file.path("filtered_barcode", "stf_4.txt"),
  stf_5 = file.path("filtered_barcode", "stf_5.txt")
) 

bclist <- lapply(bcFiles, readLines)

# Load the gene model file
gtf <- import(args$gtf)
gtf <- subset(gtf, type == "gene")
gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")

## Get expression matrices
h5Files <- c(
  stf_2 = file.path(args$data, "stf_2/outs/raw_feature_bc_matrix.h5"),
  stf_3 = file.path(args$data, "stf_3/outs/raw_feature_bc_matrix.h5"),
  stf_4 = file.path(args$data, "stf_4/outs/raw_feature_bc_matrix.h5"),
  stf_5 = file.path(args$data, "stf_5/outs/raw_feature_bc_matrix.h5")
)

# Get GEX matrices only with cells in the permissive white lists
mlist <- ReadH5mat(libs, bclist, h5Files, "Gene Expression")

# Standardize the genes available in each matrices
gene_list <- lapply(mlist, row.names) |>
  unlist() |>
  unique()
cgenes <- intersect(gene_list, gtf$gene_name)

# Remove rRNA counts
cgenes <- cgenes[!grepl("SrR", cgenes)]

mlist <- lapply(
  mlist, function(m) {
    to_pad <- setdiff(cgenes, row.names(m))
    padm <- as(matrix(0, nrow = length(to_pad), ncol = ncol(m)), "RsparseMatrix")
    row.names(padm) <- to_pad
    out <- rbind(m, padm)
    out <- out[cgenes, ]
    out <- as(out, "CsparseMatrix")
    return(out)
  }
)

m <- do.call(cbind, mlist)
rm(mlist)
saveRDS(m, "raw_count.rds")
