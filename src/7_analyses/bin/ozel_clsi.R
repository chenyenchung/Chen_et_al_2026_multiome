#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$origin_map <- "static/origin_map.csv"
  args$obj_dir <- "static/ozel_2020_obj"
  args$out_dir <- "result/ozel_2020_clsi"
  args$dims <- 120
  args$k <- 20
  args$cluster_col <- "Clustering"
}

has_arg <- function(name) {
  name %in% names(args) && !is.null(args[[name]])
}

stopifnot(
  has_arg("proot"),
  dir.exists(args[["proot"]]),
  has_arg("origin_map"),
  has_arg("out_dir")
)

if (!has_arg("clsi_helper")) {
  args[["clsi_helper"]] <- file.path(
    args[["proot"]], "src", "common", "clsi_helpers.R"
  )
}

if (has_arg("obj")) {
  args[["obj"]] <- strsplit(args[["obj"]], ",")[[1]]
}

stopifnot(
  file.exists(args[["origin_map"]]),
  file.exists(args[["clsi_helper"]]),
  !has_arg("obj_dir") || dir.exists(args[["obj_dir"]]),
  !has_arg("obj") || all(file.exists(args[["obj"]]))
)

renv::load(args[["proot"]])
source(args[["clsi_helper"]])

library(Seurat)
library(RANN)

if (!has_arg("dims")) {
  args[["dims"]] <- 120
}
if (!has_arg("k")) {
  args[["k"]] <- 20
}
if (!has_arg("cluster_col")) {
  args[["cluster_col"]] <- "Clustering"
}

dims <- as.integer(args[["dims"]])
k <- as.integer(args[["k"]])
stopifnot(
  length(dims) == 1,
  !is.na(dims),
  dims > 0,
  length(k) == 1,
  !is.na(k),
  k > 0
)

origin_map <- read.csv(args[["origin_map"]], stringsAsFactors = FALSE)
required_cols <- c("cluster_id", "tori", "nori", "sori_simplified")
missing_cols <- setdiff(required_cols, colnames(origin_map))
if (length(missing_cols) > 0) {
  stop(
    "origin map is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

origin_map$cluster_id <- as.character(origin_map$cluster_id)
origin_map$tori <- trimws(origin_map$tori)
origin_map$nori <- trimws(origin_map$nori)
origin_map$sori_simplified <- trimws(origin_map$sori_simplified)
origin_lut <- origin_map[!duplicated(origin_map$cluster_id), ]
row.names(origin_lut) <- origin_lut$cluster_id

get_object_files <- function(args) {
  if (has_arg("obj")) {
    return(args[["obj"]])
  }
  stopifnot(has_arg("obj_dir"))
  list.files(args[["obj_dir"]], pattern = "\\.rds$", full.names = TRUE)
}

process_object <- function(obj_file) {
  message("Processing ", basename(obj_file))
  obj <- readRDS(obj_file)

  meta <- obj@meta.data
  if (!args[["cluster_col"]] %in% colnames(meta)) {
    stop(
      basename(obj_file),
      " is missing metadata column ",
      args[["cluster_col"]]
    )
  }
  if (!"pca" %in% names(obj@reductions)) {
    stop(basename(obj_file), " is missing a pca reduction")
  }

  pca <- obj@reductions[["pca"]]@cell.embeddings
  if (ncol(pca) < dims) {
    stop(
      basename(obj_file),
      " has only ",
      ncol(pca),
      " PCs, but ",
      dims,
      " were requested"
    )
  }

  cluster_id <- as.character(meta[[args[["cluster_col"]]]])
  mapped <- origin_lut[cluster_id, required_cols]

  keep <- !is.na(mapped$tori) &
    mapped$tori != "" &
    !is.na(mapped$nori) &
    mapped$nori != "" &
    !is.na(mapped$sori_simplified) &
    mapped$sori_simplified != "" &
    mapped$sori_simplified != "multi"

  keep[is.na(keep)] <- FALSE
  keep_cells <- row.names(meta)[keep]
  if (length(keep_cells) < 2) {
    stop(basename(obj_file), " has fewer than two eligible cells after filtering")
  }

  pca_keep <- pca[keep_cells, seq_len(dims), drop = FALSE]
  k_use <- min(k, nrow(pca_keep))
  nn_idx <- RANN::nn2(
    data = pca_keep,
    query = pca_keep,
    k = k_use
  )$nn.idx

  out <- data.frame(
    spatial_clsi = compute_clsi(
      mapped$sori_simplified[keep],
      nn_idx,
      label_name = "sori_simplified"
    ),
    temporal_clsi = compute_clsi(
      mapped$tori[keep],
      nn_idx,
      label_name = "tori"
    ),
    notch_clsi = compute_clsi(
      mapped$nori[keep],
      nn_idx,
      label_name = "nori"
    ),
    spatial_label = mapped$sori_simplified[keep],
    temporal_label = mapped$tori[keep],
    atlas_cluster_or_cell_type = cluster_id[keep],
    row.names = keep_cells
  )

  if (!dir.exists(args[["out_dir"]])) {
    dir.create(args[["out_dir"]], recursive = TRUE)
  }
  out_file <- file.path(
    args[["out_dir"]],
    paste0(tools::file_path_sans_ext(basename(obj_file)), "_clsi.csv")
  )
  write.csv(out, out_file)
  message(
    "Wrote ",
    out_file,
    " (",
    nrow(out),
    " / ",
    nrow(meta),
    " cells)"
  )
}

obj_files <- get_object_files(args)
if (length(obj_files) == 0) {
  stop("no .rds object files found")
}

invisible(lapply(obj_files, process_object))
