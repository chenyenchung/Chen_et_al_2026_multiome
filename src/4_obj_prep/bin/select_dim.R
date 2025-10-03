#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$obj <-"int/filtered_wnn.rds"
  args$proot <- "/scratch/ycc520/thesis"
  args$ncpus <- "14"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$ncpus),
  !is.na(as.integer(args$ncpus))
)

knn_idx <- function(Z, k = 30) {
  # drop self
  nn2(Z, Z, k = k + 1)$nn.idx[, -1, drop = FALSE] 
}

jaccard_row <- function(a, b) {
  # a,b are integer vectors of length k
  length(intersect(a, b)) / length(union(a, b))
}

knn_stability <- function(Z, d1, d2, k = 30) {
  nn1 <- knn_idx(Z[, 1:d1, drop = FALSE], k = k)
  nn2 <- knn_idx(Z[, 1:d2, drop = FALSE], k = k)
  out <- sapply(seq_len(nrow(nn1)), function(i) {
    return(jaccard_row(nn1[i,], nn2[i,]))
  })
  return(out)
}

renv::load(args$proot)

library(Seurat)
library(Signac)
library(RANN)
library(future.apply)
library(tidyr)
library(pheatmap)
library(ggplot2)
library(dplyr)

plan(multisession, workers = as.integer(args$ncpus))

obj <- readRDS(args$obj)
dim_step <- 5
dim_vec <- seq(from = 50, to = 250, by = dim_step)

# GEX: Harmony-corrected PCA
hpca <- obj@reductions$hpca@cell.embeddings
gex_dim <- future_lapply(dim_vec, function(d) {
  out <- knn_stability(
    hpca,
    d1 = d - dim_step,
    d2 = d,
    k = 30
  )
  return(out)
})

gex_tbl <- do.call(cbind, gex_dim)
colnames(gex_tbl) <- paste0("d=", dim_vec)
write.csv(gex_tbl, "gex_dim_jaccard.csv", row.names = FALSE)


gex_out <- apply(gex_tbl, 2, quantile, seq(0, 1, 0.01))
row_labs <- row.names(gex_out)
row_labs <- vapply(seq_along(row_labs), function(i) {
  if (i %% 5 == 1) {
    return(row_labs[i])
  }
  return("")
}, FUN.VALUE = character(1))
pheatmap(
  gex_out, cluster_cols = FALSE, cluster_rows = FALSE,
  labels_row = row_labs,
  main = "# of PCs vs Neighborhood Similarity (k = 30)",
  filename = "PCA_jaccard_heatmap.png"
)

gex_out <- as.data.frame(gex_out)
gex_out$percentile <- row.names(gex_out)
gex_out <- pivot_longer(
  gex_out, cols = -percentile, names_to = "ndim", values_to = "jsim"
)

p <- gex_out %>%
  mutate(
    percentile = as.integer(sub("%$", "", percentile)) / 100,
    ndim = as.integer(sub("^d=", "", ndim))
  ) %>%
  filter(percentile %in% c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) %>%
  ggplot(aes(x = ndim, y = jsim)) +
  geom_smooth(
    aes(color = as.factor(percentile)), se = FALSE, method = "loess",
    formula = "y ~ x"
  ) +
  scale_x_continuous(breaks = seq(from = 50, to = 250, by = 10)) +
  scale_color_viridis_d() +
  theme_classic() +
  labs(
    x = "# of Dimensions",
    y = paste0("Jaccard Similarity\n(delta = ", dim_step, ")"),
    title = "Gene Expression (PCA)",
    color = "Percentile"
  )

ggsave("PCA_jaccard_curve.png", p)

# ATAC: Harmony-corrected LSI
hlsi <- obj@reductions$hlsi@cell.embeddings
atac_dim <- future_lapply(dim_vec, function(d) {
  out <- knn_stability(
    hlsi,
    d1 = d - dim_step,
    d2 = d,
    k = 30
  )
  return(out)
})

atac_tbl <- do.call(cbind, atac_dim)
colnames(atac_tbl) <- paste0("d=", dim_vec)
write.csv(atac_tbl, "atac_dim_jaccard.csv", row.names = FALSE)

atac_out <- apply(atac_tbl, 2, quantile, seq(0, 1, 0.01))
row_labs <- row.names(atac_out)
row_labs <- vapply(seq_along(row_labs), function(i) {
  if (i %% 5 == 1) {
    return(row_labs[i])
  }
  return("")
}, FUN.VALUE = character(1))
pheatmap(
  atac_out, cluster_cols = FALSE, cluster_rows = FALSE,
  labels_row = row_labs,
  main = "# of LSI Dim vs Neighborhood Similarity (k = 30)",
  filename = "LSI_jaccard_heatmap.png"
)

atac_out <- as.data.frame(atac_out)
atac_out$percentile <- row.names(atac_out)
atac_out <- pivot_longer(
  atac_out, cols = -percentile, names_to = "ndim", values_to = "jsim"
)
p <- atac_out %>%
  mutate(
    percentile = as.integer(sub("%$", "", percentile)) / 100,
    ndim = as.integer(sub("^d\\=", "", ndim))
  ) %>%
  filter(percentile %in% c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)) %>%
  ggplot(aes(x = ndim, y = jsim)) +
  geom_smooth(
    aes(color = as.factor(percentile)), se = FALSE, method = "loess",
    formula = "y ~ x"
  ) +
  scale_x_continuous(breaks = seq(from = 50, to = 250, by = 10)) +
  scale_color_viridis_d() +
  theme_classic() +
  labs(
    x = "# of Dimensions",
    y = paste0("Jaccard Similarity\n(delta = ", dim_step, ")"),
    title = "ATAC (LSI)",
    color = "Percentile"
  )

ggsave("LSI_jaccard_curve.png", p)
