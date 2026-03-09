#!/usr/bin/env nextflow
nextflow.preview.types = true
nextflow.enable.strict = true

params.root = '/scratch/ycc520/thesis'
params.obj = '/scratch/ycc520/thesis/int/objects/annotated.rds'
params.fasta = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
params.all_peaks = '/scratch/ycc520/thesis/int/consensus_peaks.bed'
params.motif = '/scratch/ycc520/thesis/static/CisBPDrosophilaALL.meme'
params.gtf = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf'

process RunDE {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path

  output:
  ideg: Path = file('de_identity/RNA.csv')
  idar: Path = file('de_identity/ATAC.csv')
  sdeg: Path = file('de_spatial/RNA.csv')
  sdar: Path = file('de_spatial/ATAC.csv')

  script:
  """
  run_de.R --proot ${params.root} \
    --obj ${obj}
  """
}

process GetDARbed {
  module 'r/4.5.1'
  cpus '1'
  memory '4GB'
  time '15m'

  input:
  dart: Path

  output:
  beds: Set<Path> = files("*.bed")

  script:
  """
  get_DAR_bed.R --proot ${params.root} \
    --sdar ${dart}
  """
}

process RunXSTREME {
  conda 'meme bedtools'
  cpus '1'
  time '1h'
  memory '16GB'

  input:
  bed: Path

  output:
  outdir: Set<Path> = files("${bed.baseName}/**")

  script:
  """
  bedtools getfasta -fi ${params.fasta} -bed ${bed} > peaks.fa
  bedtools getfasta -fi ${params.fasta} -bed ${params.all_peaks} > bg.fa
  xstreme --o ${bed.baseName} --p peaks.fa \
    --n bg.fa --m ${params.motif} --dna
  """
}

process Fig2ISS {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path

  output:
  sch: Set<Path> = files("fig/f2/A_*.pdf")
  iss_st: Path = file('fig/f2/clsi_spatial_vs_temporal.pdf')
  iss_st_tbl: Path = file('tbl/clsi_spatial_vs_temporal.csv')
  iss_sc: Path = file('sup_fig/s2/clsi_spatial_over_class.pdf')
  iss_sc_tbl: Path = file('tbl/clsi_spatial_over_class.csv')

  script:
  """
  fig2_ISS.R --proot ${params.root} \
    --obj ${obj}
  """
}

process Fig2DECount {
  module 'r/4.5.1'
  cpus '1'
  memory '8GB'
  time '30m'

  input:
  dar: Path
  deg: Path

  output:
  deg_st: Path = file('fig/f2/DEG_counts.pdf')
  deg_g_st: Path = file('fig/f2/DEG_counts_guide.pdf')
  dar_st: Path = file('fig/f2/DAR_counts.pdf')
  dar_g_st: Path = file('fig/f2/DAR_counts_guide.pdf')
  de_sonly: Path = file('sup_fig/s2/spatial_DE_counts.pdf')
  de_g_sonly: Path = file('sup_fig/s2/spatial_DE_counts_guide.pdf')

  script:
  """
  fig2_DE_number.R --proot ${params.root} \
    --dar ${dar} --deg ${deg}
  """
}

process Fig2DARIdeo {
  module 'r/4.5.1'
  cpus '1'
  memory '8GB'
  time '30m'

  input:
  dar: Path
  gtf: Path

  output:
  ideo: Path = file('fig/f2/dar_ideo.pdf')
  vsx: Path = file('fig/f2/dar_Vsx.pdf')
  bi: Path = file('fig/f2/dar_bi.pdf')
  optix: Path = file('sup_fig/s2/dar_optix.pdf')
  dpp: Path = file('sup_fig/s2/dar_dpp.pdf')
  rx: Path = file('sup_fig/s2/dar_rx.pdf')

  script:
  """
  fig2_differential_ideo.R --proot ${params.root} \
    --dar ${dar} --gtf ${gtf}
  """
}

process SFig2QC {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path

  output:
  comp: Path = file('sup_fig/s2/lib_comp.pdf')
  qc: Path = file('sup_fig/s2/qc_metric.pdf')

  script:
  """
  fig2_QC.R --proot ${params.root} \
    --obj ${obj}
  """
}

process SFig2UMAP {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path

  output:
  umap: Path = file('sup_fig/s2/seurat_cluster_umap.pdf')

  script:
  """
  fig2_raw_umap.R --proot ${params.root} \
    --obj ${obj}
  """
}
workflow {

  main:
  de_ch = RunDE(file(params.obj))
  bed_ch = GetDARbed(de_ch.sdar).flatten()
  xstr_ch = RunXSTREME(bed_ch)
  f2_iss_ch = Fig2ISS(file(params.obj))
  f2_de_count_ch = Fig2DECount(de_ch.idar, de_ch.ideg)
  f2_ideo_ch = Fig2DARIdeo(de_ch.idar, file(params.gtf))
  fig2_qc_ch = SFig2QC(file(params.obj))
  fig2_umap_ch = SFig2UMAP(file(params.obj))

  publish:
  id_deg = de_ch.ideg
  id_dar = de_ch.idar
  sp_deg = de_ch.sdeg
  sp_dar = de_ch.sdar
  xstreme = xstr_ch
  f2_sch = f2_iss_ch.sch
  f2_iss_st = f2_iss_ch.iss_st
  f2_iss_st_tbl = f2_iss_ch.iss_st_tbl
  f2_iss_sc = f2_iss_ch.iss_sc
  f2_iss_sc_tbl = f2_iss_ch.iss_sc_tbl
  f2_deg_st = f2_de_count_ch.deg_st
  f2_deg_st_g = f2_de_count_ch.deg_g_st
  f2_dar_st = f2_de_count_ch.dar_st
  f2_dar_st_g = f2_de_count_ch.dar_g_st
  s2_de = f2_de_count_ch.de_sonly
  s2_de_g = f2_de_count_ch.de_g_sonly
  f2_ideo = f2_ideo_ch.ideo
  f2_ideo_vsx = f2_ideo_ch.vsx
  f2_ideo_bi = f2_ideo_ch.bi
  s2_ideo_optix = f2_ideo_ch.optix
  s2_ideo_dpp = f2_ideo_ch.dpp
  s2_ideo_rx = f2_ideo_ch.rx
  s2_comp = fig2_qc_ch.comp
  s2_qc = fig2_qc_ch.qc
  s2_umap = fig2_umap_ch.umap
}

output {
  id_deg {
  }
  id_dar {
  }
  sp_deg {
  }
  sp_dar {
  }
  xstreme {
  }
  f2_sch {
  }
  f2_iss_st {
  }
  f2_iss_st_tbl {
  }
  f2_iss_sc {
  }
  f2_iss_sc_tbl {
  }
  f2_deg_st {
  }
  f2_deg_st_g {
  }
  f2_dar_st {
  }
  f2_dar_st_g {
  }
  s2_de {
  }
  s2_de_g {
  }
  f2_ideo {
  }
  f2_ideo_vsx {
  }
  f2_ideo_bi {
  }
  s2_ideo_optix {
  }
  s2_ideo_dpp {
  }
  s2_ideo_rx {
  }
  s2_comp {
  }
  s2_qc {
  }
  s2_umap {
  }
}
