#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$obj <-"int/objects/filtered_wnn.rds"
  args$resolution <- "1.1"
}

stopifnot(
  !is.null(args$obj),
  file.exists(args$obj),
  !is.null(args$proot),
  dir.exists(args$proot),
  !is.null(args$resolution),
  !is.na(as.numeric(args$resolution))
)

renv::load(args$proot)

library(Seurat)
library(Signac)
library(Matrix)
library(RSpectra)
library(bluster)
library(ggplot2)
library(igraph)
library(data.table)

obj <- readRDS(args$obj)

obj <- FindMultiModalNeighbors(
  obj, reduction.list = list("hpca", "hlsi"),
  # Selected by Jaccard Similarity Stability
  dims.list = list(1:120, 1:110)
)

# Resolution QC: Eigen Gap
k_eigs <- 50

# Get WNN shared nearest neighbor graph
single_community <- FALSE
while (!single_community) {
  A <- obj@graphs$wsnn
  Matrix::diag(A) <- 0
  
  # Disconnected component detection
  # (Where isolates were identified originally)
  g <- graph_from_adjacency_matrix(
    A, mode = "plus", weighted = TRUE, diag = FALSE
  )
  cc <- components(g)
  if (cc$no == 1) {single_community = TRUE}
  
  if (cc$no > 1) {
    message("There is more than one disjoin communties.")
    cfreq <- table(cc$membership)
    print(cfreq)
    to_trim <- names(cfreq)[
      setdiff(seq_len(cc$no), which.max(cfreq))
    ]
    message(
      "Dropping community #", paste(to_trim, collapse = ", "),
      " (size: ", cfreq[names(cfreq) %in% to_trim], ")."
    )
    bc_drop <- names(cc$membership)[as.character(cc$membership) %in% to_trim]
    obj <- subset(obj, cells = setdiff(Cells(obj), bc_drop))
    
    obj <- FindMultiModalNeighbors(
      obj, reduction.list = list("hpca", "hlsi"),
      # Selected by Jaccard Similarity Stability
      dims.list = list(1:120, 1:110)
    )
  }
}

# Resolution sweep
obj <- FindClusters(
  obj, resolution = seq(0.1, 1.2, 0.05), graph.name = "wsnn"
)

# Ensure symmetric before graph Laplacian-based QC
if (isSymmetric(A)) {
  # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
  d <- Matrix::rowSums(A)
  Dinv_sqrt <- Diagonal(x = 1 / sqrt(pmax(d, 1e-12)))
  Lsym <- Diagonal(n = nrow(A)) - Dinv_sqrt %*% A %*% Dinv_sqrt
  
  # Smallest eigenvalues of Lsym
  eig <- RSpectra::eigs_sym(Lsym, k = k_eigs + 1, which = "SM")
  lambda <- sort(Re(eig$values))
  
  # Eigengaps
  gaps <- c(diff(lambda), 0)
  
  # Plot df
  lambda_df <- data.frame(
    index = seq_along(lambda),
    lambda = lambda,
    gaps = gaps
  )
  
  p1 <- lambda_df |>
    ggplot(aes(x = index, y = lambda)) +
    geom_bar(stat = "identity", fill = NA, color = "black") +
    labs(
      y = "Eigenvalue (Normalized Laplacian)"
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank()
    )
  
  p2 <- lambda_df |>
    ggplot(aes(x = index, y = gaps)) +
    geom_bar(stat = "identity", fill = NA, color = "black") +
    labs(
      y = "Eigengap\n(lambda[k + 1]-lambda[k])"
    ) +
    geom_vline(xintercept = 36) +
    theme_classic() +
    theme(
      axis.title.x = element_blank()
    )
  ggsave("graph_lapacian_eig.pdf", p1)
  ggsave("graph_lapacian_eig_gap.pdf", p2)
  
  # Resolution QC: Silhouette
  # Dimentionality sweep
  sil_df <- lapply(
    seq(10, k_eigs, 10), function(m) {
      X <- eig$vectors[, 2:(m + 1), drop = FALSE]
      idents <- grep("^wsnn_res", colnames(obj[[]]), value = TRUE)
      sil_df <- lapply(idents, function(ident) {
        cl <- obj[[]][[ident]]
        sil <- bluster::approxSilhouette(X, clusters = cl)
        out <- data.frame(
          sil_score = sil$width,
          sil_cluster = sil$cluster,
          resolution = sub("^wsnn_res\\.", "", ident)
        )
        return(out)
      })
      sil_df <- do.call(rbind, sil_df)
      sil_df$k <- m
      return(sil_df)
    }
  )
  sil_df <- do.call(rbind, sil_df)
}

sil_dt <- as.data.table(sil_df)
p3 <- sil_dt[
  , .(avg_score = mean(sil_score), sd_score = sd(sil_score)),
  by = c("k", "resolution")
] |>
  ggplot(aes(x = as.numeric(resolution), y = avg_score)) +
  geom_pointrange(
    aes(
      color = as.factor(k),
      ymin = avg_score - sd_score,
      ymax = avg_score + sd_score
    )
  ) +
  geom_vline(xintercept = 1.1) +
  scale_color_viridis_d() +
  theme_classic() +
  labs(color = "# of Dim\n(Graph Laplacian)") +
  facet_grid(k ~ ., scales = "free_y")
ggsave("graph_lapacian_silhouette.pdf", p3, height = 12, width = 6)

Idents(obj) <- paste0("wsnn_res.", args$resolution)

# Drop unused resolutions
cluster_cols <- grep("^wsnn_res\\.", colnames(obj[[]]), value = TRUE)
for (col in cluster_cols) {
  obj@meta.data[[col]] <- NULL
}
obj[[paste0("wsnn_res.", args$resolution)]] <- Idents(obj)
obj <- RunUMAP(obj, nn.name = "weighted.nn")
saveRDS(obj, "raw_cluster.rds")
