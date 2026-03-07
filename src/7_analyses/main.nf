#!/usr/bin/env nextflow
nextflow.preview.types = true
nextflow.enable.strict = true

params.root = '/scratch/ycc520/thesis'
params.obj = '/scratch/ycc520/thesis/int/objects/annotated.rds'
params.fasta = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
params.all_peaks = '/scratch/ycc520/thesis/int/consensus_peaks.bed'
params.motif = '/scratch/ycc520/thesis/static/CisBPDrosophilaALL.meme'

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

workflow {

  main:
  de_ch = RunDE(file(params.obj))
  bed_ch = GetDARbed(de_ch.sdar).flatten()
  xstr_ch = RunXSTREME(bed_ch)
  f2_iss_ch = Fig2ISS(file(params.obj))

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
}

output {
  id_deg {}
  id_dar {}
  sp_deg {}
  sp_dar {}
  xstreme {}
  f2_sch {}
  f2_iss_st {}
  f2_iss_st_tbl {}
  f2_iss_sc {}
  f2_iss_sc_tbl {}
}
