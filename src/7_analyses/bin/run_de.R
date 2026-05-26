#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "int/objects/annotated.rds"
  args$edger_helper <- "src/7_analyses/bin/edger_ql_helpers.R"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$edger_helper),
  file.exists(args$edger_helper),
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)

library(Seurat)
library(Signac)
library(edgeR)
source(args$edger_helper)

obj <- readRDS(args$obj)


# Do neuron as a test
subset_test <- function(
  obj,
  class,
  condition,
  design_formula,
  lfc = NULL,
  mode,
  assay
) {
  # Create identities for pseudobulking
  obj$cde <- paste(obj$class, obj[[condition]][, 1], obj$library, sep = "@")

  # Get aggregated count matrix
  cmt <- AggregateExpression(
    obj,
    assay = assay,
    return.seurat = TRUE,
    group.by = "cde"
  )
  mtx <- LayerData(cmt, "counts")

  # Generate metadata from sample label
  meta <- do.call(rbind.data.frame, strsplit(cmt@meta.data$orig.ident, "@"))
  colnames(meta) <- c("class", condition, "library")
  meta$bio_rep <- "A"
  meta$bio_rep[meta$library == "stf-3"] <- "B"
  meta$bio_rep[meta$library %in% c("stf-4", "stf-5")] <- "C"
  meta[[condition]] <- factor(meta[[condition]])

  if (condition == "spatial_origin") {
    meta$spatial_origin <- factor(
      meta$spatial_origin,
      levels = c("pxb", "dpp", "optix"),
      labels = c("Vsx", "Bi", "Optix")
    )
  } else if (condition == "temporal_identity") {
    if (class != "neuron") {
      meta$temporal_identity <- factor(
        meta$temporal_identity,
        levels = c("Early", "Late")
      )
    } else {
      meta$temporal_identity <- factor(
        meta$temporal_identity,
        levels = c(
          "Hth",
          "Hth/Opa",
          "Opa/Erm",
          "Erm/Ey",
          "Ey/Hbn",
          "Hbn/Opa/Slp",
          "Slp/D",
          "D/B-H1"
        )
      )
    }
  }

  meta$bio_rep <- factor(meta$bio_rep)

  smtx <- mtx[, meta$class == class]
  smeta <- meta[meta$class == class, ]
  smeta[[condition]] <- droplevels(smeta[[condition]])
  row.names(smeta) <- colnames(smtx)

  out <- run_edger_ql(
    smtx,
    smeta,
    design_formula = design_formula,
    lfc = lfc,
    mode = mode
  )
  return(out)
}

if (!dir.exists("de_identity")) {
  dir.create("de_identity")
}

for (assay in c("RNA", "ATAC")) {
  de_tbl <- lapply(c("NE", "NB", "GMC", "neuron"), function(class) {
    spatial <- subset_test(
      obj,
      class,
      "spatial_origin",
      design_formula = ~ 0 + spatial_origin + bio_rep,
      lfc = 0.5,
      mode = "ovr",
      assay = assay
    )
    spatial$contrast_type <- "spatial"
    out <- spatial
    if (class != "NE") {
      temporal <- subset_test(
        obj,
        class,
        "temporal_identity",
        design_formula = ~ 0 + temporal_identity + bio_rep,
        lfc = 0.5,
        mode = "ovr",
        assay = assay
      )
      temporal$contrast_type <- "temporal"
      out <- rbind(out, temporal)
    }
    if (class == "neuron") {
      notch <- subset_test(
        obj,
        class,
        "notch_status",
        design_formula = ~ 0 + notch_status + bio_rep,
        lfc = 1,
        mode = "ovr",
        assay = assay
      )
      notch$contrast_type <- "notch"
      out <- rbind(out, notch)
    }
    out$subset <- class
    return(out)
  })
  de_tbl <- do.call(rbind, de_tbl)
  write.csv(
    de_tbl,
    file.path(
      "de_identity",
      paste0(assay, ".csv")
    ),
    row.names = FALSE
  )
}

if (!dir.exists("de_spatial")) {
  dir.create("de_spatial")
}

for (assay in c("RNA", "ATAC")) {
  sde_tbl <- lapply(c("NE", "NB", "GMC", "neuron"), function(class) {
    spatial <- subset_test(
      obj,
      class,
      "spatial_origin",
      design_formula = ~ 0 + spatial_origin + bio_rep,
      lfc = 0.5,
      mode = "pair",
      assay = assay
    )
    spatial$subset <- class
    return(spatial)
  })
  sde_tbl <- do.call(rbind, sde_tbl)
  write.csv(
    sde_tbl,
    file.path(
      "de_spatial",
      paste0(assay, ".csv")
    ),
    row.names = FALSE
  )
}
