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
library(ggplot2)
library(ggrastr)
library(dunn.test)
library(tidyr)

obj <- readRDS(args$obj)

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
clsi <- obj[[]][, c("class", "spatial_clsi", "temporal_clsi")]
clsi_plot <- subset(clsi, class %in% c("NB", "GMC", "neuron"))
clsi_plot <- pivot_longer(clsi_plot, cols = -class, names_to = "ident", values_to = "clsi")
clsi_plot$ident <- factor(
  clsi_plot$ident, levels = c("temporal_clsi", "spatial_clsi"),
  labels = c("Temporal", "Spatial")
)

clsi_test <- lapply(c("NB", "GMC", "neuron"), function(i) {
  df_test <- subset(clsi_plot, class == i)
  wilcox_obj <- wilcox.test(clsi ~ ident, data = df_test, conf.int = TRUE)
  out <- data.frame(
    ident = i,
    p.value = wilcox_obj$p.value,
    confint_low = wilcox_obj$conf.int[1],
    confint_up = wilcox_obj$conf.int[2]
  )
})
clsi_test <- do.call(rbind, clsi_test)
if (!dir.exists("tbl")) {dir.create("tbl")}
write.csv(
  clsi_test, file.path("tbl", "clsi_spatial_vs_temporal.csv"), row.names = FALSE
)

p_svt <- clsi_plot |>
  ggplot(aes(x = ident, y = clsi)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(aes(color = ident), outlier.size = 0.1) +
  annotate(
    geom = "segment",
    x = 1, xend = 2,
    y = 1.1, yend = 1.1
  ) +
  annotate(
    geom = "segment",
    x = 1, xend = 2,
    y = 1.2, yend = 1.2,
    color = "transparent"
  ) +
  annotate(
    geom = "text",
    x = 1.5, y = 1.1, label = "*", size = 8
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_text(margin = margin(r = 5)),
    axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  guides(color = "none") +
  scale_color_manual(values = c("#67001F", "#053061")) +
  labs(y = "Identity Segregation Score") +
  facet_grid(. ~ class)

ggsave(
  filename = file.path("fig", "f2", "clsi_spatial_vs_temporal.pdf"),
  plot = p_svt, height = 3, width = 4
)

# For supplemental figure showing ISS is increasing for spatial identity
spatial_only <- subset(clsi_plot, ident == "Spatial")
dunn_obj <- dunn.test(
  spatial_only$clsi, spatial_only$class, method = "bonferroni"
)
dunn_df <- data.frame(
  comparison = dunn_obj$comparisons,
  p.adj = dunn_obj$P.adjusted
)
write.csv(
  dunn_df, file.path("tbl", "clsi_spatial_over_class.csv"), row.names = FALSE
)

p_s_over_time <- spatial_only |>
  ggplot(aes(x = class, y = clsi)) +
  geom_hline(yintercept = 0) +
  geom_boxplot(notch = TRUE) +
  geom_smooth(aes(group = 1), method = "lm", linetype = 2, formula = y ~ x) +
  annotate(
    geom = "segment",
    x = c(1.1, 2.1, 1), xend = c(1.9, 2.9, 3),
    y = c(1.1, 1.1, 1.4), yend = c(1.1, 1.1, 1.4)
  ) +
  annotate(
    geom = "text",
    x = c(1.5, 2.5, 2), y = c(1.15, 1.15, 1.4), label = "*", size = 8
  ) +
  theme_classic() +
  theme(
    axis.title.x = element_blank()
  ) +
  labs(y = "Identity Segregation Score")

if (!dir.exists("sup_fig/s2")) {dir.create("sup_fig/s2", recursive = TRUE)}
ggsave(
  filename = file.path("sup_fig", "s2", "clsi_spatial_over_class.pdf"),
  plot = p_s_over_time, height = 3, width = 4
)
