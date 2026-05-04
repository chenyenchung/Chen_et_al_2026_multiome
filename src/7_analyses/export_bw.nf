#!/usr/bin/env nextflow
nextflow.enable.types = true
nextflow.enable.strict = true

params {
  root: String = '/scratch/ycc520/thesis'
  cluster_csv: Path = null
  bam: Path = '/scratch/ycc520/thesis/data/bam_new'
  scbamop: Path = '/scratch/ycc520/thesis/src/5_cluster_annotation/bin/scbamop'
  binsize: Integer = 1
  outdir: String = 'ATAC'
}


process PrepareSplitMeta {
  conda 'python=3.12'
  cpus '1'
  memory '2GB'
  time '15m'

  input:
  cluster_csv: Path

  output:
  meta: Set<Path> = files('meta/*.csv')

  script:
  """
  mkdir -p meta
  prepare_export_bw_meta.py \
    --input "${cluster_csv}" \
    --bam-root "${params.bam}" \
    --output-dir meta
  """
}

process SplitATACBam {
  conda 'samtools'
  cpus '1'
  memory '1GB'
  time '2h'

  input:
  tuple(lib: String, bam: Path, meta: Path)

  output:
  cbam: Set<Path> = files("${lib}/*.bam")

  script:
  """
  samtools view -H "${bam}" | \
    awk 'BEGIN{OFS="\t"} \$0 ~ /^@HD/ && !(\$0 ~ /VN:/){\$0 = sprintf("%s\t%s\t%s", \$1, "VN:1.6", \$2)}{print \$0}' > \
    header_fix.sam
  samtools reheader header_fix.sam "${bam}" > fixed.bam
  "${params.scbamop}" split -f fixed.bam -m "${meta}" -q 30 -L 18 --atac -o "${lib}"
  """
}

process GetATACbw {
  conda 'deeptools samtools'
  cpus '4'
  memory '8GB'
  time '2h'

  input:
  bam: Path

  output:
  bws: Set<Path> = files("${bam.baseName}_*.bw")

  script:
  """
  samtools index "${bam}"
  alignmentSieve -p ${task.cpus} \
    -b "${bam}" --ATACshift -o "${bam.baseName}_s.bam"

  samtools sort -@ ${task.cpus} "${bam.baseName}_s.bam" -o "${bam.baseName}_sorted.bam"
  samtools index "${bam.baseName}_sorted.bam"
  bamCoverage -b "${bam.baseName}_sorted.bam" -o "${bam.baseName}_s.bw" \
    -p ${task.cpus} --samFlagExclude 1024 \
    --Offset 1 --binSize ${params.binsize} --normalizeUsing CPM \
    --ignoreForNormalization Y dmel_mitochondrion_genome
  bamCoverage -b "${bam.baseName}_sorted.bam" -o "${bam.baseName}_e.bw" \
    -p ${task.cpus} --samFlagExclude 1024 \
    --Offset -1 --binSize ${params.binsize} --normalizeUsing CPM \
    --ignoreForNormalization Y dmel_mitochondrion_genome
  """
}

process AvgATACbw {
  conda 'deeptools=3.5.6'
  cpus '1'
  memory '2GB'
  time '30m'

  input:
  tuple(cluster: String, bws: Bag<Path>)

  stage:
  stageAs bws, 'bw??/*'

  output:
  cbw: Path = file("${cluster}_atac.bw")

  script:
  def bw_args = bws.collect { bw -> "\"${bw}\"" }.join(' ')
  """
  bigwigAverage -b ${bw_args} -o "${cluster}_atac.bw"
  """
}

workflow {
  main:
  meta_ch = PrepareSplitMeta(params.cluster_csv)
    .flatten()
    .map { meta ->
      return tuple(meta.baseName, meta)
    }

  abam_ch = meta_ch.map { lib, meta ->
    return tuple(lib, file("${params.bam}/${lib}/outs/atac_possorted_bam.bam"), meta)
  }

  cabam_ch = SplitATACBam(abam_ch)
  abw_ch = GetATACbw(cabam_ch.flatten())

  obw_ch = AvgATACbw(
    abw_ch.flatten().map { bw ->
      def cluster = bw.baseName.replaceFirst(/_(s|e)$/, '')
      return tuple(cluster, bw)
    }.groupTuple()
  )

  publish:
  atac_bw = obw_ch
}

output {
  atac_bw {
    path "bw/${params.outdir}"
  }
}
