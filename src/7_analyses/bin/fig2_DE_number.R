#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$dar <- "result/de_identity/ATAC.csv"
  args$deg <- "result/de_identity/RNA.csv"
}

stopifnot(
  !is.null(args$dar),
  file.exists(args$dar),
  !is.null(args$deg),
  file.exists(args$deg),
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)

library(ggplot2)
library(data.table)

de_bar <- function(x, type = "DEGs", by = "contrast_type") {
  dt <- copy(x)
  dt <- dt[, dir := ifelse(logFC > 0, 1, -1)][
    , dir_label := ifelse(logFC > 0, "Positive", "Negative")
  ][
    FDR < 0.05 & contrast_type != "notch", .N,
    by = c("subset", "contrast", by, "dir", "dir_label")
  ]
  
  dt$dir_label <- factor(
    dt$dir_label,
    levels = c("Positive", "Negative")
  )
  dt$subset <- factor(
    dt$subset,
    levels = c("NB", "GMC", "neuron"),
    labels = c("NB", "GMC", "Neuron")
  )
  dt$contrast <- sub(" vs rest", "", dt$contrast)
  dt$contrast <- factor(
    dt$contrast,
    levels = c("Vsx", "Optix", "Bi", "Hth", "Hth/Opa", "Opa/Erm", "Erm/Ey", 
               "Ey/Hbn", "Hbn/Opa/Slp", "Slp/D", "D/B-H1")
  )
  dt$contrast_type <- factor(
    dt$contrast_type,
    levels = c("spatial", "temporal"),
    labels = c("Spatial", "Temporal")
  )
  
  p <- dt |>
    ggplot(aes(x = .data[[by]], y = N * dir, fill = contrast)) + 
    geom_bar(stat = "identity") +
    facet_grid(dir_label ~ subset, scales = "free_y") +
    scale_fill_manual(
      values = c(
        "Vsx" = "#1B9E77",
        "Optix" = "#7570B3",
        "Bi" = "#D95F02",
        "Hth" = "#C76F6B",
        "Hth/Opa" = "#EFC8B9",
        "Opa/Erm" = "#FEE699",
        "Erm/Ey" = "#79A68C",
        "Ey/Hbn" = "#82C9C5",
        "Hbn/Opa/Slp" = "#B4C6E6",
        "Slp/D" = "#A6A1CD",
        "D/B-H1" = "#A46690"
      )
    ) +
    labs(
      y = paste0("# of ", type),
    ) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      legend.title = element_blank(),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
    )
  return(p)
}

ggsave_sep <- function(
    plot, filename,
    width = NULL, height = NULL, limitsize = FALSE
) {
  legendsp <- cowplot::get_plot_component(plot, "guide-box", return_all = TRUE)
  plot <- plot + theme(legend.position="none")
  
  if (!inherits(legendsp, "gtable")) {
    ## Drop empty elements
    to_keep <- sapply(legendsp, function(x) !"zeroGrob" %in% class(x))
    legendsp <- legendsp[to_keep]
    if (length(legendsp) == 1) {
      legendsp <- legendsp[[1]]
    } else {
      stop("Legend extraction error: There is more than 1 item.")
    }
  }
  basename <- strsplit(filename, "\\.")[[1]][1]
  ext <- strsplit(filename, "\\.")[[1]][2]
  ggsave(
    filename = filename, plot = plot,
    width = width , height = height, limitsize = limitsize
  )
  ggsave(
    filename = paste0(basename, "_guide.", ext),
    legendsp
  )
}

rna <- fread(args$deg)
rna$de_class <- "RNA"
atac <- fread(args$dar)
atac$de_class <- "ATAC"

if (!dir.exists("fig/f2")) {dir.create("fig/f2", recursive = TRUE)}
ggsave_sep(
  de_bar(rna, type = "DEGs"),
  file.path("fig", "f2", "DEG_counts.pdf"), width = 3, height = 4
)
ggsave_sep(
  de_bar(atac, type = "DARs"),
  file.path("fig", "f2", "DAR_counts.pdf"), width = 3, height = 4
)

spatial_only <- rbind(rna[contrast_type == "spatial"], atac[contrast_type == "spatial"])

if (!dir.exists("sup_fig/s2")) {dir.create("sup_fig/s2", recursive = TRUE)}
ggsave_sep(
  de_bar(spatial_only, by = "de_class", type = "DEGs/DARs"),
  file.path("sup_fig", "s2", "spatial_DE_counts.pdf"),
  width = 3, height = 4
)
