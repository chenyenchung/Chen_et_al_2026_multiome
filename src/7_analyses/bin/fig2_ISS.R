#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <- "int/objects/annotated.rds"
  args$ozel_clsi_dir <- "result/ozel_2020_clsi"
}

has_arg <- function(name) {
  name %in% names(args) && !is.null(args[[name]])
}

resolve_path <- function(path, root = args[["proot"]]) {
  if (file.exists(path) || grepl("^/", path)) {
    return(path)
  }
  return(file.path(root, path))
}

stopifnot(
  has_arg("obj"),
  has_arg("proot"),
  dir.exists(args[["proot"]])
)

args[["obj"]] <- resolve_path(args[["obj"]])
if (has_arg("ozel_clsi_dir")) {
  args[["ozel_clsi_dir"]] <- resolve_path(args[["ozel_clsi_dir"]])
}
if (has_arg("ozel_clsi")) {
  args[["ozel_clsi"]] <- vapply(
    strsplit(args[["ozel_clsi"]], ",")[[1]],
    resolve_path,
    FUN.VALUE = character(1)
  )
}

stopifnot(
  file.exists(args[["obj"]]),
  !has_arg("ozel_clsi_dir") || dir.exists(args[["ozel_clsi_dir"]]),
  !has_arg("ozel_clsi") || all(file.exists(args[["ozel_clsi"]]))
)

renv::load(args[["proot"]])

library(Seurat)
library(Signac)
library(ggplot2)
library(ggrastr)
library(dunn.test)
library(tidyr)

stage_levels <- c(
  "NE",
  "NB",
  "GMC",
  "Neuron",
  "P15",
  "P30",
  "Adult"
)
report_stage_levels <- c("NE", "NB", "GMC", "Neuron", "P15", "P30", "Adult")
plot_stage_levels <- report_stage_levels
clsi_effect_threshold <- 0.2
clsi_adjustment_method <- "Bonferroni"

get_ozel_clsi_files <- function() {
  if (has_arg("ozel_clsi")) {
    return(args[["ozel_clsi"]])
  }
  if (has_arg("ozel_clsi_dir")) {
    return(list.files(
      args[["ozel_clsi_dir"]],
      pattern = "_clsi\\.csv$",
      full.names = TRUE
    ))
  }
  character(0)
}

read_ozel_clsi <- function(files) {
  stage_order <- c("P15", "P30", "Adult")
  files_by_stage <- setNames(rep(NA_character_, length(stage_order)), stage_order)
  file_stage <- sub("_clsi\\.csv$", "", basename(files))
  files_by_stage[file_stage] <- files

  missing_stage <- names(files_by_stage)[is.na(files_by_stage)]
  if (length(missing_stage) > 0) {
    stop(
      "missing Ozel cLSI CSVs for stages: ",
      paste(missing_stage, collapse = ", ")
    )
  }

  ozel_clsi <- lapply(stage_order, function(stage) {
    clsi <- read.csv(
      files_by_stage[[stage]], row.names = 1, stringsAsFactors = FALSE
    )
    clsi$cell_id <- row.names(clsi)
    clsi$dataset <- "Ozel_2020"
    clsi$stage <- stage
    clsi$object_class <- "neuron"
    clsi$library_or_sample <- stage
    optional_cols <- c(
      "spatial_label", "temporal_label", "atlas_cluster_or_cell_type"
    )
    for (column in setdiff(optional_cols, colnames(clsi))) {
      clsi[[column]] <- NA_character_
    }
    clsi
  })
  do.call(rbind, ozel_clsi)
}

format_report_p <- function(p_value) {
  if (is.na(p_value)) {
    return(NA_character_)
  }
  if (p_value == 0) {
    return(paste0("< ", format(.Machine$double.xmin, scientific = TRUE)))
  }
  format(p_value, scientific = TRUE, digits = 6, trim = TRUE)
}

obj <- readRDS(args[["obj"]])

# UMAPs for the workflow schematic
plotdf <- Embeddings(obj, reduction = "umap") |>
  as.data.frame()

plotdf <- merge(plotdf, obj[["spatial_origin"]], by = "row.names")
for (i in unique(plotdf$spatial_origin)) {
  p <- subset(plotdf, spatial_origin == i) |>
    ggplot(aes(x = umap_1, y = umap_2, color = spatial_origin)) +
    geom_point_rast(size = 0.2) +
    guides(color = "none") +
    scale_color_manual(
      values = c(
        "pxb" = "#1B9E77",
        "optix" = "#7570B3",
        "dpp" = "#D95F02"
      ),
      na.value = "grey80",
      breaks = c("pxb", "optix", "dpp")
    ) +
    theme_void()
  if (!dir.exists("fig/f2")) {dir.create("fig/f2", recursive = TRUE)}
  ggsave(
    filename = file.path("fig", "f2", paste0("A_", i, ".pdf")),
    plot = p, width = 5, height = 4
  )
}

# Identity segregation score comparison and raw-data export
l3_required_cols <- c(
  "class", "library", "spatial_clsi", "temporal_clsi",
  "spatial_origin", "temporal_identity", "nn_pred_label"
)
l3_missing_cols <- setdiff(l3_required_cols, colnames(obj[[]]))
if (length(l3_missing_cols) > 0) {
  stop(
    "annotated object is missing required metadata columns: ",
    paste(l3_missing_cols, collapse = ", ")
  )
}
l3_clsi <- obj[[]][, l3_required_cols]
l3_clsi <- subset(l3_clsi, class %in% c("NE", "NB", "GMC", "neuron"))
l3_clsi$cell_id <- row.names(l3_clsi)
l3_clsi$dataset <- "L3_multiome"
l3_clsi$stage <- as.character(l3_clsi$class)
l3_clsi$stage[l3_clsi$stage == "neuron"] <- "Neuron"
l3_clsi$object_class <- as.character(l3_clsi$class)
l3_clsi$library_or_sample <- as.character(l3_clsi$library)
l3_clsi$spatial_label <- as.character(l3_clsi$spatial_origin)
l3_clsi$temporal_label <- as.character(l3_clsi$temporal_identity)
l3_clsi$atlas_cluster_or_cell_type <- as.character(l3_clsi$nn_pred_label)
l3_clsi <- l3_clsi[, c(
  "cell_id", "dataset", "stage", "object_class", "library_or_sample",
  "spatial_clsi", "temporal_clsi", "spatial_label", "temporal_label",
  "atlas_cluster_or_cell_type"
)]

ozel_files <- get_ozel_clsi_files()
if (length(ozel_files) == 0) {
  stop("no Ozel cLSI CSVs supplied")
}
ozel_clsi <- read_ozel_clsi(ozel_files)
ozel_clsi <- ozel_clsi[, colnames(l3_clsi)]

clsi <- rbind(l3_clsi, ozel_clsi)
clsi$stage <- factor(clsi$stage, levels = stage_levels)
clsi <- subset(clsi, !is.na(stage))

raw_iss <- subset(clsi, stage %in% report_stage_levels)
raw_iss$stage <- as.character(raw_iss$stage)
raw_iss <- raw_iss[, c(
  "cell_id", "dataset", "stage", "object_class",
  "spatial_clsi", "temporal_clsi", "spatial_label", "temporal_label",
  "atlas_cluster_or_cell_type"
)]
colnames(raw_iss)[colnames(raw_iss) == "spatial_clsi"] <- "spatial_ISS"
colnames(raw_iss)[colnames(raw_iss) == "temporal_clsi"] <- "temporal_ISS"
if (!dir.exists("tbl")) {dir.create("tbl")}
write.csv(raw_iss, file.path("tbl", "ISS_raw.csv"), row.names = FALSE, na = "")

clsi_plot <- pivot_longer(
  clsi,
  cols = c("temporal_clsi", "spatial_clsi"),
  names_to = "ident",
  values_to = "clsi"
)
clsi_plot <- subset(clsi_plot, !is.na(clsi))
clsi_plot$ident <- factor(
  clsi_plot$ident, levels = c("temporal_clsi", "spatial_clsi"),
  labels = c("Temporal", "Spatial")
)
clsi_plot$stage <- factor(clsi_plot$stage, levels = stage_levels)

clsi_test <- lapply(report_stage_levels, function(i) {
  df_test <- subset(clsi, stage == i)
  n_cells <- sum(!is.na(df_test$spatial_clsi))
  n_libraries_or_samples <- length(unique(
    df_test$library_or_sample[!is.na(df_test$spatial_clsi)]
  ))
  median_spatial <- if (n_cells > 0) {
    median(df_test$spatial_clsi, na.rm = TRUE)
  } else {
    NA_real_
  }
  df_paired <- subset(
    df_test, !is.na(spatial_clsi) & !is.na(temporal_clsi)
  )
  if (i == "NE") {
    return(data.frame(
      stage = i,
      n_cells = n_cells,
      n_libraries_or_samples = n_libraries_or_samples,
      median_spatial = median_spatial,
      median_temporal = NA_real_,
      median_paired_difference = NA_real_,
      HL_location_shift = NA_real_,
      HL_CI_low = NA_real_,
      HL_CI_high = NA_real_,
      minimum_effect_threshold = NA_real_,
      p_raw_numeric = NA_real_,
      interpretation = paste(
        "Spatial ISS only; temporal ISS and paired comparison not performed"
      )
    ))
  }
  if (nrow(df_paired) == 0) {
    stop("no paired spatial and temporal ISS values for stage ", i)
  }
  clsi_diff <- df_paired$spatial_clsi - df_paired$temporal_clsi
  wilcox_obj <- wilcox.test(
    clsi_diff, mu = 0, conf.int = TRUE, exact = FALSE
  )
  wilcox_greater <- wilcox.test(
    clsi_diff, mu = clsi_effect_threshold,
    alternative = "greater", exact = FALSE
  )
  wilcox_less <- wilcox.test(
    clsi_diff, mu = -clsi_effect_threshold,
    alternative = "less", exact = FALSE
  )
  p_effect <- min(
    1,
    2 * min(wilcox_greater$p.value, wilcox_less$p.value)
  )
  out <- data.frame(
    stage = i,
    n_cells = nrow(df_paired),
    n_libraries_or_samples = length(unique(df_paired$library_or_sample)),
    median_spatial = median(df_paired$spatial_clsi),
    median_temporal = median(df_paired$temporal_clsi),
    median_paired_difference = median(clsi_diff),
    HL_location_shift = unname(wilcox_obj$estimate),
    HL_CI_low = unname(wilcox_obj$conf.int[1]),
    HL_CI_high = unname(wilcox_obj$conf.int[2]),
    minimum_effect_threshold = clsi_effect_threshold,
    p_raw_numeric = p_effect,
    interpretation = NA_character_
  )
})
clsi_test <- do.call(rbind, clsi_test)
test_idx <- !is.na(clsi_test$p_raw_numeric)
clsi_test$p_adjusted_numeric <- NA_real_
clsi_test$p_adjusted_numeric[test_idx] <- p.adjust(
  clsi_test$p_raw_numeric[test_idx], method = "bonferroni"
)
clsi_test$signif <- NA_character_
clsi_test$signif[test_idx] <- ifelse(
  clsi_test$p_adjusted_numeric[test_idx] < 0.001, "***",
  ifelse(
    clsi_test$p_adjusted_numeric[test_idx] < 0.01, "**",
    ifelse(clsi_test$p_adjusted_numeric[test_idx] < 0.05, "*", "n.s.")
  )
)
clsi_test$interpretation[test_idx] <- ifelse(
  clsi_test$p_adjusted_numeric[test_idx] < 0.05,
  ifelse(
    clsi_test$HL_location_shift[test_idx] > 0,
    "Spatial ISS exceeds temporal ISS by at least the minimum effect threshold",
    "Temporal ISS exceeds spatial ISS by at least the minimum effect threshold"
  ),
  "No evidence that the spatial-temporal ISS shift exceeds the minimum effect threshold"
)
clsi_test$p_raw <- vapply(
  clsi_test$p_raw_numeric, format_report_p, FUN.VALUE = character(1)
)
clsi_test$p_adjusted <- vapply(
  clsi_test$p_adjusted_numeric, format_report_p, FUN.VALUE = character(1)
)
clsi_test$adjustment_method <- ifelse(
  test_idx, clsi_adjustment_method, "not_performed"
)
clsi_report <- clsi_test[, c(
  "stage", "n_cells", "n_libraries_or_samples", "median_spatial",
  "median_temporal", "median_paired_difference", "HL_location_shift",
  "HL_CI_low", "HL_CI_high", "minimum_effect_threshold", "p_raw",
  "p_adjusted", "adjustment_method", "interpretation"
)]
write.csv(
  clsi_report,
  file.path("tbl", "clsi_spatial_vs_temporal.csv"),
  row.names = FALSE,
  na = ""
)

clsi_plot_visible <- subset(clsi_plot, stage %in% plot_stage_levels)
clsi_plot_visible$stage <- factor(
  clsi_plot_visible$stage, levels = plot_stage_levels
)

clsi_sig <- subset(
  clsi_test,
  stage %in% plot_stage_levels & stage != "NE"
)
clsi_sig$stage <- factor(clsi_sig$stage, levels = plot_stage_levels)
clsi_y <- aggregate(clsi ~ stage, clsi_plot_visible, max)
clsi_sig <- merge(clsi_sig, clsi_y, by = "stage", all.x = TRUE)
clsi_range <- diff(range(clsi_plot_visible$clsi, na.rm = TRUE))
if (clsi_range == 0) {clsi_range <- 1}
clsi_sig$label_y <- clsi_sig$clsi + clsi_range * 0.08

p_svt <- clsi_plot_visible |>
  ggplot(aes(x = stage, y = clsi, color = ident, fill = ident)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(
    outlier.size = 0.1,
    position = position_dodge2(width = 0.75, preserve = "single"),
    alpha = 0.4
  ) +
  geom_text(
    data = clsi_sig,
    aes(x = stage, y = label_y, label = signif),
    inherit.aes = FALSE,
    size = 3
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    legend.title = element_blank()
  ) +
  scale_color_manual(values = c("Temporal" = "#67001F", "Spatial" = "#053061")) +
  scale_fill_manual(values = c("Temporal" = "#67001F", "Spatial" = "#053061")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  guides(fill = "none", color = "none") +
  labs(y = "Identity Segregation Score")

ggsave(
  filename = file.path("fig", "f2", "clsi_spatial_vs_temporal.pdf"),
  plot = p_svt, height = 3.5, width = 7.5
)

# For supplemental figure showing ISS is increasing for spatial identity
spatial_only <- subset(clsi_plot, ident == "Spatial")
spatial_only_visible <- subset(spatial_only, stage %in% plot_stage_levels)
spatial_only_visible$stage <- factor(
  spatial_only_visible$stage, levels = plot_stage_levels
)
dunn_obj <- dunn.test(
  spatial_only$clsi, spatial_only$stage, method = "bonferroni"
)
dunn_df <- data.frame(
  comparison = dunn_obj$comparisons,
  p.adj = dunn_obj$P.adjusted
)
write.csv(
  dunn_df, file.path("tbl", "clsi_spatial_over_class.csv"), row.names = FALSE
)

p_s_over_time <- spatial_only_visible |>
  ggplot(aes(x = stage, y = clsi)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(notch = TRUE) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1)
  ) +
  labs(y = "Identity Segregation Score")

if (!dir.exists("sup_fig/s2")) {dir.create("sup_fig/s2", recursive = TRUE)}
ggsave(
  filename = file.path("sup_fig", "s2", "clsi_spatial_over_class.pdf"),
  plot = p_s_over_time, height = 3.5, width = 5
)
