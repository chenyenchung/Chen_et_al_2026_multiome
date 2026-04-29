make_primary_col_map <- function(meta, design_formula, primary_factor) {
  if (anyNA(meta[[primary_factor]])) {
    stop("primary_factor contains NA values: ", primary_factor)
  }

  if (is.factor(meta[[primary_factor]])) {
    meta[[primary_factor]] <- droplevels(meta[[primary_factor]])
  } else {
    meta[[primary_factor]] <- factor(
      meta[[primary_factor]],
      levels = unique(meta[[primary_factor]])
    )
  }

  levs <- levels(meta[[primary_factor]])
  if (length(levs) < 2) {
    stop("primary_factor has <2 levels; can't do contrasts.")
  }

  design <- model.matrix(design_formula, data = meta)
  mm_primary <- model.matrix(
    reformulate(primary_factor, intercept = FALSE),
    data = meta
  )

  if (!identical(rownames(design), rownames(meta))) {
    stop("Design matrix rows do not match metadata rows. Check for NA values in design covariates.")
  }
  if (!identical(rownames(mm_primary), rownames(meta))) {
    stop("Primary-factor model matrix rows do not match metadata rows.")
  }
  if (anyDuplicated(colnames(design))) {
    stop("Design matrix has duplicated column names; contrasts would be ambiguous.")
  }
  if (anyDuplicated(colnames(mm_primary))) {
    stop("Primary-factor model matrix has duplicated column names; contrasts would be ambiguous.")
  }

  primary_cols <- setNames(character(length(levs)), levs)
  for (lvl in levs) {
    expected <- as.integer(meta[[primary_factor]] == lvl)
    matches <- vapply(
      colnames(mm_primary),
      function(col) {
        observed <- as.integer(mm_primary[, col] != 0)
        identical(observed, expected)
      },
      logical(1)
    )

    if (sum(matches) != 1) {
      stop(
        "Could not uniquely map primary-factor level to design column: ",
        lvl
      )
    }

    primary_cols[[lvl]] <- colnames(mm_primary)[matches]
  }

  missing_cols <- setdiff(unname(primary_cols), colnames(design))
  if (length(missing_cols) > 0) {
    stop(
      "Primary-factor columns are not all found in design: ",
      paste(missing_cols, collapse = ", "),
      ". Use a no-intercept design for the primary factor, e.g. ~ 0 + ",
      primary_factor,
      " + covariate."
    )
  }
  for (lvl in names(primary_cols)) {
    expected <- as.integer(meta[[primary_factor]] == lvl)
    observed <- as.integer(design[, primary_cols[[lvl]]] != 0)

    if (!identical(observed, expected)) {
      stop(
        "Primary-factor design column does not encode expected level: ",
        lvl
      )
    }
  }

  list(meta = meta, design = design, primary_cols = primary_cols)
}

check_contrast_inputs <- function(design, primary_cols) {
  if (
    is.null(names(primary_cols)) ||
      anyNA(names(primary_cols)) ||
      any(names(primary_cols) == "")
  ) {
    stop("primary_cols must be named by primary-factor level.")
  }
  if (anyDuplicated(names(primary_cols))) {
    stop("primary_cols has duplicated level names.")
  }
  if (anyDuplicated(unname(primary_cols))) {
    stop("primary_cols maps multiple levels to the same design column.")
  }
  if (anyDuplicated(colnames(design))) {
    stop("Design matrix has duplicated column names; contrasts would be ambiguous.")
  }
  missing_cols <- setdiff(unname(primary_cols), colnames(design))
  if (length(missing_cols) > 0) {
    stop(
      "Contrast columns are not present in design: ",
      paste(missing_cols, collapse = ", ")
    )
  }
}

make_ovr_contrast <- function(focal_level, design, primary_cols) {
  check_contrast_inputs(design, primary_cols)

  levs <- names(primary_cols)
  cn <- colnames(design)
  v <- setNames(rep(0, length(cn)), cn)

  focal_col <- primary_cols[[focal_level]]
  if (is.null(focal_col)) {
    stop("No design column found for focal level: ", focal_level)
  }
  v[focal_col] <- 1

  other_levs <- setdiff(levs, focal_level)
  if (length(other_levs) == 0) {
    stop("No 'rest' levels found for ", focal_level)
  }

  other_cols <- unname(primary_cols[other_levs])
  if (anyNA(other_cols)) {
    stop("No design column found for one or more rest levels.")
  }

  v[other_cols] <- -(1 / length(other_levs))
  return(v)
}

make_pairwise_contrast <- function(design, primary_cols) {
  check_contrast_inputs(design, primary_cols)

  levs <- names(primary_cols)
  cn <- colnames(design)
  lev_comb <- combn(levs, 2)

  out <- list()
  for (i in seq_len(ncol(lev_comb))) {
    lvl1 <- lev_comb[1, i]
    lvl2 <- lev_comb[2, i]
    col1 <- primary_cols[[lvl1]]
    col2 <- primary_cols[[lvl2]]

    if (is.null(col1) || is.null(col2)) {
      stop("No design column found for contrast: ", lvl1, " vs ", lvl2)
    }

    v <- setNames(rep(0, length(cn)), cn)
    v[col1] <- 1
    v[col2] <- -1
    out[[paste(lvl1, "vs", lvl2)]] <- v
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

  primary_map <- make_primary_col_map(meta, design_formula, primary_factor)
  meta <- primary_map$meta
  design <- primary_map$design
  primary_cols <- primary_map$primary_cols
  levs <- names(primary_cols)

  # edgeR workflow
  y <- DGEList(counts = counts)
  y <- calcNormFactors(y, method = "TMM")

  keep <- filterByExpr(
    y,
    design = design,
    min.count = min.count,
    min.total.count = min.total.count
  )
  y <- y[keep, , keep.lib.sizes = FALSE]

  y <- estimateDisp(y, design, robust = robust)
  fit <- glmQLFit(y, design, robust = robust, prior.count = prior_count)

  # Build one-vs-rest contrasts for every level of primary_factor
  if (mode == "ovr") {
    contrast_list <- setNames(
      lapply(levs, make_ovr_contrast, design, primary_cols),
      levs
    )
  } else {
    contrast_list <- make_pairwise_contrast(design, primary_cols)
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
  attr(out, "primary_cols") <- primary_cols
  row.names(out) <- NULL
  return(out)
}
