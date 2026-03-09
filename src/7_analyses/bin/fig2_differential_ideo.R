#!/usr/bin/env Rscript
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

if (interactive()) {
  args$proot <- "/scratch/ycc520/thesis"
  args$dar <- "result/de_identity/ATAC.csv"
  args$gtf <- "/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf"
}


stopifnot(
  !is.null(args$dar),
  file.exists(args$dar),
  !is.null(args$gtf),
  file.exists(args$gtf),
  !is.null(args$proot),
  dir.exists(args$proot)
)

renv::load(args$proot)

library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
library(GenomeInfoDb)
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)

plot_mkideo <- function(
    spa_gr, tem_gr, ext = 1e5, nlayers = 10, tstep = 0.3, wsize = 5e4,
    saveas = NULL, class_pal = NULL, mode = "density", genes = FALSE,
    gene_size = 0.5, chr = "auto", zoom = NULL, lsize = 1, klwd = 0.01,
    width = 5, height = 4
) {
  stopifnot(mode %in% c("density", "region"))
  if (is.null(class_pal)) {
    class_pal <- c(
      "NB" = "#377EB8",
      "GMC" = "#4DAF4A",
      "neuron" = "#984EA3"
    )
  }
  
  .plot <- function(
    spa_gr, tem_gr, ext, nlayers, tstep, wsize, class_pal, mode, genes,
    gene_size, chr, zoom, klwd, lsize
  ) {
    pp <- getDefaultPlotParams(plot.type = 3)
    if (genes) {
      pp$data2height <- pp$data1height / 2
    }
    kp <- plotKaryotype(
      genome = "dm6", plot.type = 3,
      chromosomes = chr, zoom = zoom, cex = lsize, lwd = klwd, plot.params = pp
    )
    kpAddChromosomeSeparators(kp, lwd = 2)
    base <- 0
    
    for (i in c("NB", "GMC", "neuron")) {
      if (mode == "region") {
        kpPlotRegions(
          kp, data = extendRegions(subset(spa_gr, class == i), extend.end = ext),
          r0 = base, r1 = base + tstep, num.layers = nlayers,
          col = class_pal[i],
          data.panel = 1
        )
      }
      if (mode == "density") {
        kpPlotDensity(
          kp, data = subset(spa_gr, class == i),
          r0 = base, r1 = base + tstep, window.size = wsize,
          col = class_pal[i],
          data.panel = 1
        )
      }
      
      kpAddLabels(
        kp, labels = i, srt = 90, pos = 1, label.margin = 0.04,
        r0 = base, r1 = base + 0.3, data.panel = 1, cex = lsize
      )
      
      if (genes) {
        if (i == "neuron") {
          genes.data <- makeGenesDataFromTxDb(
            txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene, karyoplot = kp
          )
          kpPlotGenes(
            kp, data = genes.data, data.panel = 2, r0 = 0, r1 = gene_size,
            add.gene.names = FALSE
          )
        }
      } else {
        if (mode == "region") {
          kpPlotRegions(
            kp, data = extendRegions(subset(tem_gr, class == i), extend.end = ext),
            r0 = base, r1 = base + tstep, num.layers = nlayers,
            col = class_pal[i],
            data.panel = 2
          )
        }
        
        if (mode == "density") {
          kpPlotDensity(
            kp, data = subset(tem_gr, class == i), window.size = wsize,
            r0 = base, r1 = base + tstep,
            col = class_pal[i],
            data.panel = 2
          )
        }
        
        kpAddLabels(
          kp, labels = i, srt = 90, pos = 1, label.margin = 0.04, cex = lsize,
          r0 = base, r1 = base + 0.27, data.panel = 2
        )
      }
      base <- base + 0.3
    }
    
  }

  if (!is.null(saveas)) {
    pdf(saveas, width = width, height = height)
    .plot(
      spa_gr, tem_gr, ext, nlayers, tstep, wsize, class_pal, mode, genes,
      gene_size, chr, zoom, klwd, lsize
    )
    dev.off()
  } else {
    .plot(
      spa_gr, tem_gr, ext, nlayers, tstep, wsize, class_pal, mode, genes,
      gene_size, chr, zoom, klwd, lsize
    )
  }
}

gtf <- import(args$gtf)
gtf <- subset(gtf, type == "gene")
seqlevelsStyle(gtf) <- "UCSC"


atac <- read.csv(args$dar)

spa <- subset(atac, contrast_type == "spatial" & FDR < 0.05)
tem <- subset(atac, contrast_type == "temporal" & FDR < 0.05)

spa_gr <- spa$feature |>
  sub("-", ":", x = _) |>
  paste0("chr", x = _) |>
  GRanges()
spa_gr$pval <- spa$PValue
spa_gr$class <- spa$subset

tem_gr <- tem$feature |>
  sub("-", ":", x = _) |>
  paste0("chr", x = _) |>
  GRanges()
tem_gr$pval <- tem$PValue
tem_gr$class <- tem$subset

if (!dir.exists("fig/f2")) dir.create("fig/f2", recursive = TRUE)
if (!dir.exists("sup_fig/s2")) dir.create("sup_fig/s2", recursive = TRUE)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.25, mode = "density",
  saveas = file.path("fig", "f2", "dar_ideo.pdf"),
  width = 6, height = 6
)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.3, nlayers = 6,
  zoom = GRanges("chrX:5510000-5630000"), mode = "region", ext = 0,
  genes = TRUE, klwd = 1, lsize = 0.5,
  saveas = file.path("fig", "f2", "dar_Vsx.pdf"),
  width = 2, height = 4.5
)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.3, nlayers = 6,
  zoom = GRanges("chrX:4350000-4550000"), mode = "region", ext = 0,
  genes = TRUE, klwd = 1, lsize = 0.5,
  saveas = file.path("fig", "f2", "dar_bi.pdf"),
  width = 2, height = 4.5
)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.3, nlayers = 6,
  zoom = GRanges("chr2R:8020000-8050000"), mode = "region", ext = 0,
  genes = TRUE, klwd = 1, lsize = 0.5,
  saveas = file.path("sup_fig", "s2", "dar_optix.pdf"),
  width = 2, height = 3
)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.3, nlayers = 6,
  zoom = GRanges("chr2L:2420000-2470000"), mode = "region", ext = 0,
  genes = TRUE, klwd = 1, lsize = 0.5,
  saveas = file.path("sup_fig", "s2", "dar_dpp.pdf"),
  width = 2, height = 3
)

plot_mkideo(
  spa_gr, tem_gr, tstep = 0.3, nlayers = 6,
  zoom = GRanges("chr2R:20910000-20940000"), mode = "region", ext = 0,
  genes = TRUE, klwd = 1, lsize = 0.5, gene_size = 0.2,
  saveas = file.path("sup_fig", "s2", "dar_rx.pdf"),
  width = 2, height = 3
)
