#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "int/objects/annotated.rds"
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
library(edgeR)

obj <- readRDS(args$obj)

make_ovr_contrast <- function(focal_level, design, primary_cols, levs) {
  cn <- colnames(design)
  v <- setNames(rep(0, length(cn)), cn)

  # focal coef is the column for that level in the no-intercept primary matrix
  focal_col <- primary_cols[levs == focal_level]
  v[focal_col] <- 1

  other_levs <- setdiff(levs, focal_level)
  other_cols <- primary_cols[match(other_levs, levs)]

  if (length(other_cols) == 0) {
    stop("No 'rest' levels found for ", focal_level)
  }

  w <- rep(1 / length(other_levs), length(other_levs)) # equal-weight rest

  v[other_cols] <- -w
  return(v)
}

make_pairwise_contrast <- function(design, primary_cols, levs) {
  cn <- colnames(design)

  lev_comb <- combn(levs, 2)

  out <- list()
  for (i in seq_len(ncol(lev_comb))) {
    v <- setNames(rep(0, length(cn)), cn)
    col1 <- primary_cols[levs == lev_comb[1, i]]
    col2 <- primary_cols[levs == lev_comb[2, i]]
    v[col1] <- 1
    v[col2] <- -1
    out[[paste(lev_comb[1, i], "vs", lev_comb[2, i])]] <- v
  }
  return(out)
}

run_edger_ql <- function(
  counts,
  meta,
  design_formula = ~ 0 + spatial_origin + batch,
  lfc = 0,
  min.count = 10,
  min.total.count = 15,
  prior_count = 1,
  primary_factor = NULL,
  mode = "pair",
  robust = TRUE
) {
  stopifnot(all(colnames(counts) == rownames(meta)))
  stopifnot(mode %in% c("pair", "ovr"))

  # Infer primary factor (first term in formula that is a column in meta)
  if (is.null(primary_factor)) {
    tl <- attr(terms(design_formula), "term.labels")
    primary_factor <- tl[tl %in% colnames(meta)][1]
    if (is.na(primary_factor) || is.null(primary_factor)) {
      stop("Couldn't infer primary_factor. Pass primary_factor explicitly.")
    }
  }
  if (!primary_factor %in% colnames(meta)) {
    stop("primary_factor not found in meta: ", primary_factor)
  }

  # Build design matrix
  design <- model.matrix(design_formula, data = meta)

  # Sanity: make sure we can find the primary factor columns in the design
  mm_primary <- model.matrix(
    reformulate(primary_factor, intercept = FALSE),
    data = meta
  )
  primary_cols <- colnames(mm_primary)

  if (!all(primary_cols %in% colnames(design))) {
    stop(
      "Primary-factor columns are not all found in design."
    )
  }

  # edgeR workflow
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMM")

  keep <- filterByExpr(
    y,
    design = design,
    min.count = min.count,
    min.total.count =
  )
  y <- y[keep, , keep.lib.sizes = FALSE]

  y <- estimateDisp(y, design, robust = robust)
  fit <- glmQLFit(y, design, robust = robust, prior.count = prior_count)

  # Build one-vs-rest contrasts for every level of primary_factor
  if (is.factor(meta[[primary_factor]])) {
    levs <- levels(meta[[primary_factor]])
  } else {
    levs <- unique(meta[[primary_factor]])
  }
  if (length(levs) < 2) {
    stop("primary_factor has <2 levels; can't do contrasts.")
  }

  if (mode == "ovr") {
    contrast_list <- setNames(
      lapply(levs, make_ovr_contrast, design, primary_cols, levs),
      levs
    )
  } else {
    contrast_list <- make_pairwise_contrast(design, primary_cols, levs)
  }

  # Test per contrast
  out_list <- lapply(names(contrast_list), function(lvl) {
    cvec <- contrast_list[[lvl]]
    if (lfc > 0) {
      res <- glmTreat(fit, contrast = cvec, lfc = lfc)
    } else {
      res <- glmQLFTest(fit, contrast = cvec)
    }
    tab <- topTags(res, n = Inf)$table
    if (mode == "ovr") {
      tab$contrast <- paste0(lvl, " vs rest")
    } else {
      tab$contrast <- lvl
    }
    tab$feature <- rownames(tab)
    return(tab)
  })

  out <- do.call(rbind, out_list)
  attr(out, "design") <- design
  attr(out, "contrasts") <- contrast_list
  row.names(out) <- NULL
  return(out)
}


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
  de_tbl <- lapply(c("NB", "GMC", "neuron"), function(class) {
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
    out <- rbind(spatial, temporal)
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
