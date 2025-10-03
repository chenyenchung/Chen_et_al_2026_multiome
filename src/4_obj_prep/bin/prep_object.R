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

ReadPeakMat <- function(path) {
  # 1. Read only the necessary columns (4 through 8)
  # We assume the file has no header
  dt <- fread(path, header = FALSE, select = c(4, 5, 6, 7, 8), 
              col.names = c("cell_name", "count", "chr", "start", "end"))
  
  # 2. Construct the row names (Feature ID)
  # Converting BED (0-based) to 1-based coordinates
  dt[, feature_id := paste0(chr, ":", start + 1, "-", end)]
  
  # 3. Convert identifiers to factors
  # This creates integer indices required for the sparse matrix
  dt[, cell_idx := factor(cell_name)]
  dt[, feature_idx := factor(feature_id)]
  
  # 4. Construct the sparse matrix
  # sparseMatrix automatically sums 'x' values for duplicate (i, j) pairs
  sp_mat <- sparseMatrix(
    i = as.integer(dt$feature_idx),
    j = as.integer(dt$cell_idx),
    x = dt$count,
    dimnames = list(levels(dt$feature_idx), levels(dt$cell_idx))
  )
  return(sp_mat)
}


library(R.utils)
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$cbc <- "int/permissive_cbc/"
  args$gtf <- "/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf"
  args$data <- "/scratch/ycc520/thesis/data/bam_new"
  args$demux <- "/scratch/ycc520/thesis/int/demux"
  args$freemux <- "/scratch/ycc520/thesis/int/freemux"
  args$out <- "base_obj.rds"
}

stopifnot(
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$cbc),
  dir.exists(args$cbc),
  !is.null(args$data),
  dir.exists(args$data),
  !is.null(args$demux),
  dir.exists(args$demux),
  !is.null(args$freemux),
  dir.exists(args$freemux),
  !is.null(args$gtf),
  file.exists(args$gtf)
)

renv::load(args$proot)
library(Seurat)
library(Matrix)
library(Signac)
library(rtracklayer)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Dmelanogaster.BDGP.dm6)
library(data.table)
library(hdf5r)

libs <- c("stf_2", "stf_3", "stf_4", "stf_5")
names(libs) <- libs

## Load permissive barcode lists
bcFiles <- c(
  stf_2 = file.path(args$cbc, "stf_2.txt"),
  stf_3 = file.path(args$cbc, "stf_3.txt"),
  stf_4 = file.path(args$cbc, "stf_4.txt"),
  stf_5 = file.path(args$cbc, "stf_5.txt")
) 

bclist <- lapply(bcFiles, readLines)

# Load the gene model file
gtf <- import(args$gtf)
gtf <- subset(gtf, type == "gene")
gtf <- keepStandardChromosomes(gtf, pruning.mode = "coarse")

## Get expression matrices
h5Files <- c(
  stf_2 = file.path(args$data, "stf_2/outs/gex_cellbender.h5"),
  stf_3 = file.path(args$data, "stf_3/outs/gex_cellbender.h5"),
  stf_4 = file.path(args$data, "stf_4/outs/gex_cellbender.h5"),
  stf_5 = file.path(args$data, "stf_5/outs/gex_cellbender.h5")
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

insCountFiles <- list(
  "stf_2" = "stf_2.tsv.gz",
  "stf_3" = "stf_3.tsv.gz",
  "stf_4" = "stf_4.tsv.gz",
  "stf_5" = "stf_5.tsv.gz"
)

plist <- lapply(libs, function(lib) {
  mat <- ReadPeakMat(insCountFiles[[lib]])
  colnames(mat) <- paste(lib, colnames(mat), sep = "#")
  return(mat)
})

# Standardize the peaks available in each matrices
peak_list <- lapply(plist, row.names) |>
  unlist() |>
  unique() |>
  GRanges() |>
  sort() |> 
  GRangesToString(grange = _, sep = c(":", "-"))

plist <- lapply(
  plist, function(m) {
    m <- as(m, "RsparseMatrix")
    to_pad <- setdiff(peak_list, row.names(m))
    padm <- as(matrix(0, nrow = length(to_pad), ncol = ncol(m)), "RsparseMatrix")
    row.names(padm) <- to_pad
    out <- rbind(m, padm)
    out <- out[peak_list, ]
    out <- as(out, "CsparseMatrix")
    return(out)
  }
)

p <- do.call(cbind, plist)
rm(plist)

# Ensure both assays have the same cells
cells_use <- intersect(colnames(m), colnames(p))
m <- m[, cells_use]
p <- p[, cells_use]

# Get demux calls
demux_path <- args$demux
meta <- list.files(demux_path)
names(meta) <- sub("\\.best", "", meta)

demux_type <- lapply(libs, function(lib) {
  calls <- read.delim(file.path(demux_path, meta[lib]))
  row.names(calls) <- calls$BARCODE
  
  # Rename to match object
  calls <- calls[bclist[[lib]], ]
  out <- calls[["DROPLET.TYPE"]]
  names(out) <- paste(lib, calls[["BARCODE"]], sep = "#")
  return(out)
})

demux_call <- lapply(libs, function(lib) {
  calls <- read.delim(file.path(demux_path, meta[lib]))
  row.names(calls) <- calls$BARCODE
  
  # Rename to match object
  calls <- calls[bclist[[lib]], ]
  out <- calls[["SNG.BEST.GUESS"]]
  names(out) <- paste(lib, calls[["BARCODE"]], sep = "#")
  return(out)
})

# Get freemux calls
freemux_path <- args$freemux
fmeta <- list.files(freemux_path)
names(fmeta) <- sub("\\..*$", "", fmeta)

freemux_type <- lapply(libs, function(lib) {
  calls <- read.delim(file.path(freemux_path, fmeta[lib]))
  row.names(calls) <- calls$BARCODE
  
  # Reorder to match ArchR object
  calls <- calls[bclist[[lib]], ]
  out <- calls[["DROPLET.TYPE"]]
  names(out) <- paste(lib, calls[["BARCODE"]], sep = "#")
  return(out)
})

freemux_call <- lapply(libs, function(lib) {
  calls <- read.delim(file.path(freemux_path, fmeta[lib]))
  row.names(calls) <- calls$BARCODE
  
  # Reorder to match ArchR object
  calls <- calls[bclist[[lib]], ]
  out <- calls[["SNG.BEST.GUESS"]]
  names(out) <- paste(lib, calls[["BARCODE"]], sep = "#")
  return(out)
})

# Stringent: Considered singlet in both methods
# Permissive: Considered singlet in one method
stringent_call <- lapply(
  libs, function(lib) {
    dtype <- demux_type[[lib]] == "SNG"
    ftype <- freemux_type[[lib]] == "SNG"
    return(dtype & ftype)
  }
)

permissive_call <- lapply(
  libs, function(lib) {
    dtype <- demux_type[[lib]] == "SNG"
    ftype <- freemux_type[[lib]] == "SNG"
    return(dtype | ftype)
  }
)

## Also consider calling consistency:
## We assume that the majority of the calls will match.
## i.e., most of the freemux label X will match demux label X'
## Under this assumption, besides droplet type calls, we will
## also examine the consistency of labels and only keep
## *the droplets that maintain label consistency between the two methods*.
consistent_calls <- lapply(
  libs, function(lib) {
    x <- demux_call[[lib]]
    y <- freemux_call[[lib]]
    stopifnot(all(names(x) == names(y)))
    tbl <- table(x, y)
    fcall <- rownames(tbl)[apply(tbl, 2, which.max)]
    names(fcall) <- colnames(tbl)
    dcall <- colnames(tbl)[apply(tbl, 1, which.max)]
    names(dcall) <- row.names(tbl)
    fcon <- fcall[as.character(y)] == x
    dcon <- dcall[as.character(x)] == y
    out <- fcon & dcon
    names(out) <- names(x)
    return(out)
  }
)

permission_call_con <- lapply(
  libs, function(lib) {
    x <- consistent_calls[[lib]]
    y <- permissive_call[[lib]]
    stopifnot(all(names(x) == names(y)))
    return(x & y)
  }
)

stringent_call_con <- lapply(
  libs, function(lib) {
    x <- consistent_calls[[lib]]
    y <- stringent_call[[lib]]
    stopifnot(all(names(x) == names(y)))
    return(x & y)
  }
)

permission_call_con_v <- Reduce(c, permission_call_con)[cells_use]
stringent_call_con_v <- Reduce(c, stringent_call_con)[cells_use]

meta_label <- Reduce(c, demux_call)[cells_use]
stopifnot(all(!is.na(meta_label)))

dgrp_line_v <- vapply(
  strsplit(meta_label, "@"), function(x) x[[1]],
  FUN.VALUE = character(1)
)
stopifnot(all(names(dgrp_line_v) == cells_use))

sori_line_v <- vapply(
  strsplit(meta_label, "@"), function(x) x[[2]],
  FUN.VALUE = character(1)
)
stopifnot(all(names(sori_line_v) == cells_use))


## Add the following to the metadata:
## 1. Permissive filter
## 2. Stringent filter
## 3. DGRP lines
## 4. Spatial origin

if (is.null(args$meta)) {
  obj <- CreateSeuratObject(
    counts = m,
    assay = "RNA",
    project = "stf"
  )
} else {
  stopifnot(file.exists(args$meta))
  meta <- read.csv(args$meta)
  obj <- CreateSeuratObject(
    counts = m,
    assay = "RNA",
    project = "stf",
    meta.data = meta
  )
}

rm(m)

obj[["ATAC"]] <- CreateChromatinAssay(
  counts = p,
  genome = seqinfo(BSgenome.Dmelanogaster.BDGP.dm6),
  sep = c(":", "-")
)
rm(p)

if (is.null(args$meta)) {
  DefaultAssay(obj) <- "RNA"
  obj <- AddMetaData(obj, permission_call_con_v, col.name = "permissive_filter")
  obj <- AddMetaData(obj, stringent_call_con_v, col.name = "stringent_filter")
  obj <- AddMetaData(obj, dgrp_line_v, col.name = "dgrp_line")
  obj <- AddMetaData(obj, sori_line_v, col.name = "spatial_origin")
  obj <- AddMetaData(obj, sori_line_v, col.name = "spatial_origin")
  obj <- subset(obj, nCount_RNA > 200L & nCount_ATAC > 300L)
  obj$library <- vapply(
    strsplit(row.names(obj@meta.data), "#"),
    function(x) {return(x[[1]])}, FUN.VALUE = character(1)
  )
}


saveRDS(obj, args$out)
