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
  "Neuron (L3)",
  "Neuron (P15)",
  "Neuron (P30)",
  "Neuron (P40)",
  "Neuron (P50)",
  "Neuron (P70)",
  "Neuron (Adult)"
)
plot_stage_levels <- setdiff(
  stage_levels,
  c("Neuron (P30)", "Neuron (P40)", "Neuron (P70)")
)
clsi_effect_threshold <- 0.2

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
  stage_order <- c("P15", "P30", "P40", "P50", "P70", "Adult")
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
    clsi <- read.csv(files_by_stage[[stage]], row.names = 1)
    clsi$stage <- paste0("Neuron (", stage, ")")
    clsi
  })
  do.call(rbind, ozel_clsi)
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

# Identity segregation score comparison
l3_clsi <- obj[[]][, c("class", "spatial_clsi", "temporal_clsi")]
l3_clsi <- subset(l3_clsi, class %in% c("NE", "NB", "GMC", "neuron"))
l3_clsi$stage <- as.character(l3_clsi$class)
l3_clsi$stage[l3_clsi$stage == "neuron"] <- "Neuron (L3)"
l3_clsi$class <- NULL

ozel_files <- get_ozel_clsi_files()
if (length(ozel_files) == 0) {
  stop("no Ozel cLSI CSVs supplied")
}
ozel_clsi <- read_ozel_clsi(ozel_files)
ozel_clsi <- ozel_clsi[, c("stage", "spatial_clsi", "temporal_clsi")]

clsi <- rbind(l3_clsi, ozel_clsi)
clsi$stage <- factor(clsi$stage, levels = stage_levels)
clsi <- subset(clsi, !is.na(stage))

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

clsi_test <- lapply(stage_levels, function(i) {
  df_test <- subset(clsi, stage == i)
  df_test <- subset(
    df_test,
    !is.na(spatial_clsi) & !is.na(temporal_clsi)
  )
  if (i == "NE" || nrow(df_test) == 0) {
    return(data.frame(
      stage = i,
      p.value = NA_real_,
      confint_low = NA_real_,
      confint_up = NA_real_,
      median_spatial = NA_real_,
      median_temporal = NA_real_,
      median_diff = NA_real_
    ))
  }
  clsi_diff <- df_test$spatial_clsi - df_test$temporal_clsi
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
    p.value = p_effect,
    confint_low = wilcox_obj$conf.int[1],
    confint_up = wilcox_obj$conf.int[2],
    median_spatial = median(df_test$spatial_clsi),
    median_temporal = median(df_test$temporal_clsi),
    median_diff = median(clsi_diff)
  )
})
clsi_test <- do.call(rbind, clsi_test)
test_idx <- !is.na(clsi_test$p.value)
clsi_test$p.adj <- NA_real_
clsi_test$p.adj[test_idx] <- p.adjust(
  clsi_test$p.value[test_idx], method = "bonferroni"
)
clsi_test$signif <- NA_character_
clsi_test$signif[test_idx] <- ifelse(
  clsi_test$p.adj[test_idx] < 0.001, "***",
  ifelse(
    clsi_test$p.adj[test_idx] < 0.01, "**",
    ifelse(clsi_test$p.adj[test_idx] < 0.05, "*", "n.s.")
  )
)
if (!dir.exists("tbl")) {dir.create("tbl")}
write.csv(
  clsi_test, file.path("tbl", "clsi_spatial_vs_temporal.csv"), row.names = FALSE
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
