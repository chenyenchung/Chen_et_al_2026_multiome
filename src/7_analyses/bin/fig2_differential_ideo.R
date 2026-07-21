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

add_gene_strand_marks <- function(
  karyoplot,
  regs,
  transcript.height,
  mid.transcript,
  mark.height,
  mark.width,
  mark.distance,
  marks.col,
  lwd,
  r0,
  r1,
  ymin,
  ymax,
  data.panel,
  clipping
) {
  if (length(regs) == 0 || mark.height == 0) {
    return(invisible(NULL))
  }

  st <- as.character(strand(regs)[1])
  if (clipping) {
    regs <- intersect(regs, karyoplot$plot.region, ignore.strand = TRUE)
    strand(regs) <- st
  }
  if (length(regs) == 0 || st == "*") {
    return(invisible(NULL))
  }

  karyoplot$beginKpPlot()
  on.exit(karyoplot$endKpPlot())

  chr <- as.character(seqnames(regs))[1]
  pp <- karyoploteR:::prepareParameters2(
    karyoplot,
    function.name = "add_gene_strand_marks",
    data = NULL,
    chr = rep(chr, 3),
    x = c(0, 0, 0),
    y = c(mid.transcript, 0, transcript.height),
    r0 = r0,
    r1 = r1,
    ymin = ymin,
    ymax = ymax,
    data.panel = data.panel
  )
  mid.transcript <- pp$y[1]
  transcript.height <- abs(pp$y[2] - pp$y[3])

  ccf <- karyoplot$coord.change.function
  mark_height <- transcript.height * mark.height / 2
  plot_mark_height <- abs(diff(ccf(
    chr = c(chr, chr),
    x = c(0, 0),
    y = c(0, mark_height),
    data.panel = data.panel
  )$y))

  asp_usr <- diff(graphics::par("usr")[c(1, 2)]) /
    diff(graphics::par("usr")[c(3, 4)])
  asp_pin <- graphics::par("pin")[2] / graphics::par("pin")[1]
  plot_mark_width <- plot_mark_height * asp_usr * asp_pin * mark.width

  for (nr in seq_len(length(regs))) {
    reg <- regs[nr]
    plot_coords <- ccf(
      chr = c(chr, chr),
      x = c(start(reg), end(reg)),
      y = c(mid.transcript, mid.transcript),
      data.panel = data.panel
    )
    reg_width <- diff(plot_coords$x)
    if (reg_width <= 1.4 * plot_mark_width) {
      next
    }

    num_marks <- (reg_width - plot_mark_width) %/%
      (plot_mark_width * mark.distance) + 1
    left_margin <- (
      reg_width - ((num_marks - 1) * plot_mark_width * mark.distance) -
        plot_mark_width
    ) / 2
    mark_starts <- plot_coords$x[1] + left_margin + plot_mark_width / 2 +
      plot_mark_width * mark.distance * seq.int(0, num_marks - 1)

    if (st == "-") {
      graphics::segments(
        x0 = mark_starts - plot_mark_width / 2,
        x1 = mark_starts + plot_mark_width / 2,
        y0 = plot_coords$y[1],
        y1 = plot_coords$y[1] + plot_mark_height,
        col = marks.col,
        lwd = lwd
      )
      graphics::segments(
        x0 = mark_starts - plot_mark_width / 2,
        x1 = mark_starts + plot_mark_width / 2,
        y0 = plot_coords$y[1],
        y1 = plot_coords$y[1] - plot_mark_height,
        col = marks.col,
        lwd = lwd
      )
    } else {
      graphics::segments(
        x0 = mark_starts - plot_mark_width / 2,
        x1 = mark_starts + plot_mark_width / 2,
        y0 = plot_coords$y[1] + plot_mark_height,
        y1 = plot_coords$y[1],
        col = marks.col,
        lwd = lwd
      )
      graphics::segments(
        x0 = mark_starts - plot_mark_width / 2,
        x1 = mark_starts + plot_mark_width / 2,
        y0 = plot_coords$y[1] - plot_mark_height,
        y1 = plot_coords$y[1],
        col = marks.col,
        lwd = lwd
      )
    }
  }

  invisible(NULL)
}

plot_thin_gene_models <- function(
  kp,
  genes.data,
  data.panel,
  r0,
  r1,
  gene_lwd,
  gene_mark_height,
  gene_mark_width,
  gene_mark_distance,
  gene_col,
  clipping = TRUE,
  non.coding.exons.height = 0.5
) {
  kp <- kpPlotGenes(
    kp,
    data = genes.data,
    data.panel = data.panel,
    r0 = r0,
    r1 = r1,
    add.gene.names = FALSE,
    add.transcript.names = FALSE,
    add.strand.marks = FALSE,
    col = NA,
    border = NA,
    gene.col = NA,
    gene.border.col = NA,
    transcript.col = NA,
    transcript.border.col = NA,
    coding.exons.col = NA,
    coding.exons.border.col = NA,
    non.coding.exons.col = NA,
    non.coding.exons.border.col = NA,
    introns.col = NA,
    marks.col = NA
  )

  transcript_pos <- kp$latest.plot$computed.values$transcript.vertical.position
  transcripts <- unlist(GRangesList(genes.data$transcripts))
  if (length(transcripts) == 0) {
    return(invisible(kp))
  }
  names(transcripts) <- as.character(transcripts$tx_id)

  for (tx_id in names(transcripts)) {
    pos <- transcript_pos[[tx_id]]
    if (is.null(pos)) {
      next
    }

    transcript <- transcripts[tx_id]
    t_y0 <- unname(pos$y0)
    t_y1 <- unname(pos$y1)
    mid_transcript <- t_y0 + (t_y1 - t_y0) / 2
    transcript_height <- t_y1 - t_y0

    coding_exons <- genes.data$coding.exons[[tx_id]]
    non_coding_exons <- genes.data$non.coding.exons[[tx_id]]
    if (is.null(coding_exons)) {
      coding_exons <- GRanges()
    }
    if (is.null(non_coding_exons)) {
      non_coding_exons <- GRanges()
    }

    all_exons <- c(coding_exons, non_coding_exons, ignore.mcols = TRUE)
    introns <- setdiff(transcript, all_exons, ignore.strand = TRUE)
    strand(introns) <- strand(transcript)

    kpSegments(
      kp,
      data = introns,
      y0 = mid_transcript,
      y1 = mid_transcript,
      col = gene_col,
      ymin = 0,
      ymax = 1,
      r0 = r0,
      r1 = r1,
      data.panel = data.panel,
      clipping = clipping,
      lwd = gene_lwd
    )
    add_gene_strand_marks(
      kp,
      introns,
      transcript.height = transcript_height,
      mid.transcript = mid_transcript,
      mark.height = gene_mark_height,
      mark.width = gene_mark_width,
      mark.distance = gene_mark_distance,
      marks.col = gene_col,
      lwd = gene_lwd,
      r0 = r0,
      r1 = r1,
      ymin = 0,
      ymax = 1,
      data.panel = data.panel,
      clipping = clipping
    )

    nc_y0 <- mid_transcript -
      (t_y1 - t_y0) * non.coding.exons.height / 2
    nc_y1 <- mid_transcript +
      (t_y1 - t_y0) * non.coding.exons.height / 2
    kpRect(
      kp,
      data = non_coding_exons,
      y0 = nc_y0,
      y1 = nc_y1,
      col = gene_col,
      border = gene_col,
      ymin = 0,
      ymax = 1,
      r0 = r0,
      r1 = r1,
      data.panel = data.panel,
      clipping = clipping,
      lwd = gene_lwd
    )
    kpRect(
      kp,
      data = coding_exons,
      y0 = t_y0,
      y1 = t_y1,
      col = gene_col,
      border = gene_col,
      ymin = 0,
      ymax = 1,
      r0 = r0,
      r1 = r1,
      data.panel = data.panel,
      clipping = clipping,
      lwd = gene_lwd
    )
  }

  invisible(kp)
}

plot_mkideo <- function(
  spa_gr,
  tem_gr,
  ext = 1e5,
  nlayers = 10,
  tstep = 0.22,
  wsize = 5e4,
  saveas = NULL,
  class_pal = NULL,
  mode = "density",
  genes = FALSE,
  gene_size = 0.5,
  gene_lwd = 0.25,
  gene_mark_height = 0.12,
  gene_mark_width = 1,
  gene_mark_distance = 16,
  gene_col = "black",
  chr = "auto",
  zoom = NULL,
  lsize = 1,
  klwd = 0.01,
  width = 5,
  height = 4,
  density_scale = c("global", "panel"),
  markers = NULL,
  marker_panel = 1,
  marker_r0 = 0,
  marker_r1 = 1,
  marker_y = 0.75
) {
  stopifnot(mode %in% c("density", "region"))
  stopifnot(is.null(markers) || methods::is(markers, "GRanges"))
  density_scale <- match.arg(density_scale)
  if (is.null(class_pal)) {
    class_pal <- c(
      "NE" = "#E41A1C",
      "NB" = "#377EB8",
      "GMC" = "#4DAF4A",
      "neuron" = "#984EA3"
    )
  }
  get_marker_labels <- function(markers) {
    label_col <- intersect(c("labels", "label"), names(mcols(markers)))
    if (length(label_col) > 0) {
      return(as.character(mcols(markers)[[label_col[1]]]))
    }
    marker_names <- names(markers)
    if (!is.null(marker_names)) {
      return(marker_names)
    }
    stop("markers must have a 'labels' metadata column or names.")
  }

  .plot <- function(
    spa_gr,
    tem_gr,
    ext,
    nlayers,
    tstep,
    wsize,
    class_pal,
    mode,
    genes,
    gene_size,
    gene_lwd,
    gene_mark_height,
    gene_mark_width,
    gene_mark_distance,
    gene_col,
    chr,
    zoom,
    klwd,
    lsize,
    density_scale,
    markers,
    marker_panel,
    marker_r0,
    marker_r1,
    marker_y
  ) {
    pp <- getDefaultPlotParams(plot.type = 3)
    if (genes) {
      pp$data2height <- pp$data1height / 2
    }
    kp <- plotKaryotype(
      genome = "dm6",
      plot.type = 3,
      chromosomes = chr,
      zoom = zoom,
      cex = lsize,
      lwd = klwd,
      plot.params = pp
    )
    kpAddChromosomeSeparators(kp, lwd = 2)
    base <- 0
    classes <- c("NE", "NB", "GMC", "neuron")
    class_step <- 1 / length(classes)
    track_height <- min(tstep, class_step)

    if (mode == "density") {
      windows <- tileGenome(
        kp$chromosome.length,
        tilewidth = wsize,
        cut.last.tile.in.chrom = TRUE
      )
      max_density <- function(gr) {
        if (length(gr) == 0) {
          return(0)
        }
        max(countOverlaps(windows, gr))
      }

      density_ymax <- c(
        spa = max(vapply(
          classes,
          function(i) max_density(subset(spa_gr, class == i)),
          numeric(1)
        )),
        tem = max(vapply(
          classes,
          function(i) max_density(subset(tem_gr, class == i)),
          numeric(1)
        ))
      )
      if (density_scale == "global") {
        density_ymax[] <- max(density_ymax)
      }
      density_ymax[density_ymax == 0] <- 1
      add_density_axis <- function(data.panel, r0, r1, ymax) {
        ymax <- unname(ymax)
        kpAxis(
          kp,
          ymin = 0,
          ymax = ymax,
          r0 = r0,
          r1 = r1,
          data.panel = data.panel,
          side = 2,
          tick.pos = c(0, ymax),
          labels = format(c(0, ymax), trim = TRUE, scientific = FALSE),
          col = "gray30",
          text.col = "gray30",
          cex = lsize * 0.5
        )
      }
    }

    for (i in classes) {
      spa_i <- subset(spa_gr, class == i)
      tem_i <- subset(tem_gr, class == i)
      if (mode == "region" && length(spa_i) > 0) {
        kpPlotRegions(
          kp,
          data = extendRegions(spa_i, extend.end = ext),
          r0 = base,
          r1 = base + track_height,
          num.layers = nlayers,
          col = class_pal[i],
          data.panel = 1
        )
      }
      if (mode == "density" && length(spa_i) > 0) {
        kpPlotDensity(
          kp,
          data = spa_i,
          r0 = base,
          r1 = base + track_height,
          window.size = wsize,
          col = class_pal[i],
          data.panel = 1,
          ymin = 0,
          ymax = density_ymax["spa"]
        )
        add_density_axis(
          data.panel = 1,
          r0 = base,
          r1 = base + track_height,
          ymax = density_ymax["spa"]
        )
      }

      kpAddLabels(
        kp,
        labels = i,
        srt = 90,
        pos = 1,
        label.margin = 0.04,
        r0 = base,
        r1 = base + class_step,
        data.panel = 1,
        cex = lsize
      )

      if (genes) {
        if (i == "neuron") {
          genes.data <- makeGenesDataFromTxDb(
            txdb = TxDb.Dmelanogaster.UCSC.dm6.ensGene,
            karyoplot = kp
          )
          kp <- plot_thin_gene_models(
            kp,
            genes.data = genes.data,
            data.panel = 2,
            r0 = 0,
            r1 = gene_size,
            gene_lwd = gene_lwd,
            gene_mark_height = gene_mark_height,
            gene_mark_width = gene_mark_width,
            gene_mark_distance = gene_mark_distance,
            gene_col = gene_col
          )
        }
      } else {
        if (mode == "region" && length(tem_i) > 0) {
          kpPlotRegions(
            kp,
            data = extendRegions(tem_i, extend.end = ext),
            r0 = base,
            r1 = base + track_height,
            num.layers = nlayers,
            col = class_pal[i],
            data.panel = 2
          )
        }

        if (mode == "density" && length(tem_i) > 0) {
          kpPlotDensity(
            kp,
            data = tem_i,
            window.size = wsize,
            r0 = base,
            r1 = base + track_height,
            col = class_pal[i],
            data.panel = 2,
            ymin = 0,
            ymax = density_ymax["tem"]
          )
          add_density_axis(
            data.panel = 2,
            r0 = base,
            r1 = base + track_height,
            ymax = density_ymax["tem"]
          )
        }

        kpAddLabels(
          kp,
          labels = i,
          srt = 90,
          pos = 1,
          label.margin = 0.04,
          cex = lsize,
          r0 = base,
          r1 = base + class_step,
          data.panel = 2
        )
      }
      base <- base + class_step
    }

    if (!is.null(markers) && length(markers) > 0) {
      kpPlotMarkers(
        kp,
        data = markers,
        labels = get_marker_labels(markers),
        data.panel = marker_panel,
        r0 = marker_r0,
        r1 = marker_r1,
        y = marker_y,
        cex = lsize * 0.7
      )
    }
  }

  if (!is.null(saveas)) {
    pdf(saveas, width = width, height = height)
    .plot(
      spa_gr,
      tem_gr,
      ext,
      nlayers,
      tstep,
      wsize,
      class_pal,
      mode,
      genes,
      gene_size,
      gene_lwd,
      gene_mark_height,
      gene_mark_width,
      gene_mark_distance,
      gene_col,
      chr,
      zoom,
      klwd,
      lsize,
      density_scale,
      markers,
      marker_panel,
      marker_r0,
      marker_r1,
      marker_y
    )
    dev.off()
  } else {
    .plot(
      spa_gr,
      tem_gr,
      ext,
      nlayers,
      tstep,
      wsize,
      class_pal,
      mode,
      genes,
      gene_size,
      gene_lwd,
      gene_mark_height,
      gene_mark_width,
      gene_mark_distance,
      gene_col,
      chr,
      zoom,
      klwd,
      lsize,
      density_scale,
      markers,
      marker_panel,
      marker_r0,
      marker_r1,
      marker_y
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

spa_gr <- unique(spa_gr)
tem_gr <- unique(tem_gr)

if (!dir.exists("fig/f2")) {
  dir.create("fig/f2", recursive = TRUE)
}
if (!dir.exists("sup_fig/s2")) {
  dir.create("sup_fig/s2", recursive = TRUE)
}

regions_of_interest <- GRanges(c(
  "chrX:5510000-5630000",
  "chrX:4350000-4550000",
  "chr2R:8020000-8050000",
  "chr2L:2420000-2470000",
  "chr2R:20910000-20940000",
  "chrX: 7307939-7311132"
))
regions_of_interest$label <- c(
  "Vsx",
  "Bi",
  "Optix",
  "Dpp",
  "Rx",
  "Brk"
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  mode = "density",
  markers = regions_of_interest,
  saveas = file.path("fig", "f2", "dar_ideo.pdf"),
  width = 6,
  height = 6
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chrX:5510000-5630000"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  saveas = file.path("fig", "f2", "dar_Vsx.pdf"),
  width = 2,
  height = 4.5
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chrX:4350000-4550000"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  saveas = file.path("fig", "f2", "dar_bi.pdf"),
  width = 2,
  height = 4.5
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chr2R:8020000-8050000"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  saveas = file.path("sup_fig", "s2", "dar_optix.pdf"),
  width = 4,
  height = 3
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chr2L:2420000-2470000"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  saveas = file.path("sup_fig", "s2", "dar_dpp.pdf"),
  width = 4,
  height = 3
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chr2R:20910000-20940000"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  gene_size = 0.2,
  saveas = file.path("sup_fig", "s2", "dar_rx.pdf"),
  width = 4,
  height = 3
)

plot_mkideo(
  spa_gr,
  tem_gr,
  tstep = 0.22,
  nlayers = 3,
  zoom = GRanges("chrX: 7307939-7311132"),
  mode = "region",
  ext = 0,
  genes = TRUE,
  klwd = 1,
  lsize = 0.5,
  gene_size = 0.7,
  saveas = file.path("sup_fig", "s2", "dar_brk.pdf"),
  width = 4,
  height = 3
)
