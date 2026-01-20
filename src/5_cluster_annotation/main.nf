#!/usr/bin/env nextflow
nextflow.preview.types = true
nextflow.enable.strict = true

params.root = '/scratch/ycc520/thesis'
params.initobj = '/scratch/ycc520/thesis/int/objects/raw_cluster.rds'
params.frags = '/scratch/ycc520/thesis/int/insertion_np'
params.bam = '/scratch/ycc520/thesis/data/bam_new'
params.gtf = '/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf'
params.fai = '/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.fai'
params.cbc = '/scratch/ycc520/thesis/int/filtered_barcode/'
params.binsize = 1 

process GetDepth {
  conda 'python=3.14'
  cpus '1'
  memory '2GB'
  time '30m'

  input:
  _frags: List<Path>

  output:
  insertion_depth: Path = file('ins_depth.csv')

  script:
  """
  for i in insertion_np/*.tsv.gz; do
    ins_depth.py -f \${i} >> 'ins_depth.csv'
  done
  """
}

process Annotate {
  module 'r/gcc/4.5.0'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path
  depth: Path
  _frags: List<Path>

  output:
  wobj: Path = file('obj/annotated.rds')
  onobj: Path = file('obj/non_obj.rds')
  offobj: Path = file('obj/noff_obj.rds')
  metadata: Path = file('metadata.csv')
  figs: Set<Path> = files('fig/**')
  supfigs: Set<Path> = files('sup/**')

  script:
  """
  annotate.R --proot ${params.root} \
    --obj ${obj} \
    --depth ${depth}
  """
}

process GetRawCount {
  module 'r/gcc/4.5.0'
  cpus '1'
  memory '16GB'
  time '1h'

  input:
  data: List<Path>
  gtf: Path
  _cbc: List<Path>

  output:
  mtx: Path = file('raw_count.rds')

  script:
  """
  get_raw.R --proot ${params.root} \
    --data ${data} \
    --gtf ${gtf}
  """
}

process GetSplitMeta {
  module 'r/gcc/4.5.0'
  cpus '1'
  memory '4GB'
  time '15m'

  input:
  (lib, split, metadata): Tuple<String, String, Path>

  output:
  Path = file("${lib}.csv")

  script:
  """
  get_metadata.R --proot ${params.root} \
    --file ${metadata} \
    --lib ${lib} \
    --split ${split} \
    --debug
  """
}

process SplitATACBam {
  module 'samtools/intel/1.20'
  cpus '1'
  memory '1GB'
  time '2h'

  input:
  (lib, bam, meta): Tuple<String, Path, Path>

  output:
  cbam: Set<Path> = files("${lib}/*.bam")

  script:
  """
  samtools view -H ${bam} |\
    awk 'BEGIN{OFS="\t"} \$0 ~ /^@HD/ && !(\$0 ~ /VN:/){\$0 = sprintf("%s\t%s\t%s", \$1, "VN:1.6", \$2)}{print \$0}' >\
    header_fix.sam
  samtools reheader header_fix.sam ${bam} > fixed.bam
  scbamop split -f fixed.bam -m ${meta} -q 30 -L 18 --atac -o ${lib}
  """
}

process SplitGexBam {
  cpus '1'
  memory '10GB'
  time '4h'

  input:
  (lib, bam, meta): Tuple<String, Path, Path>

  output:
  cbam: Set<Path> = files("${lib}/*.bam")

  script:
  """
  scbamop split -f ${bam} -m ${meta} -L 18 -l 12 -o ${lib} -d -q 30
  """
}

process GetATACbw {
  module 'deeptools/3.5.0:samtools/intel/1.20'
  cpus '4'
  memory '8GB'
  time '2h'

  input:
  bam: Path

  output:
  bws: Set<Path> = files("${bam.baseName}_*.bw")

  script:
  """
  samtools index ${bam}
  alignmentSieve -p ${task.cpus} \
    -b ${bam} --ATACshift -o ${bam.baseName}_s.bam

  samtools sort -@ ${task.cpus} ${bam.baseName}_s.bam -o ${bam.baseName}_sorted.bam
  samtools index ${bam.baseName}_sorted.bam
  bamCoverage -b ${bam.baseName}_sorted.bam -o ${bam.baseName}_s.bw \
    -p ${task.cpus} \
    --Offset 1 --binSize ${params.binsize} --normalizeUsing CPM \
    --ignoreForNormalization Y dmel_mitochondrion_genome
  bamCoverage -b ${bam.baseName}_sorted.bam -o ${bam.baseName}_e.bw \
    -p ${task.cpus} \
    --Offset -1 --binSize ${params.binsize} --normalizeUsing CPM \
    --ignoreForNormalization Y dmel_mitochondrion_genome
  """
}

process GetGEXbw {
  module 'deeptools/3.5.0:samtools/intel/1.20'
  cpus '1'
  memory '8GB'
  time '1h'

  input:
  bam: Path

  output:
  bw: Path = file("${bam.baseName}.bw")

  script:
  """
  samtools index ${bam}
  bamCoverage -b ${bam} -o ${bam.baseName}.bw \
    --binSize ${params.binsize} --normalizeUsing CPM \
    --ignoreForNormalization Y dmel_mitochondrion_genome
  """
}

process AvgGEXbw {
  conda 'deeptools==3.5.6'
  cpus '1'
  memory '2GB'
  time '30m'

  input:
  (cluster, bws): Tuple<String, Bag<Path>>

  stage:
  stageAs 'input.bw', bws

  output:
  cbw: Path = file("${cluster}_gex.bw")

  script:
  """
  bigwigAverage -b ${bws} -o ${cluster}_gex.bw
  """
}

process AvgATACbw {
  conda 'deeptools==3.5.6'
  cpus '1'
  memory '2GB'
  time '30m'

  input:
  (cluster, bws): Tuple<String, Bag<Path>>

  stage:
  stageAs 'input.bw', bws

  output:
  cbw: Path = file("${cluster}_atac.bw")

  script:
  """
  bigwigAverage -b ${bws} -o ${cluster}_atac.bw
  """
}

process FeatureSelection {
  module 'r/gcc/4.5.0'
  cpus '1'
  memory '32GB'
  time '4h'

  input:
  obj: Path

  output:
  g: Path = file('varg.rds')
  p: Path = file('varp.rds')

  script:
  """
  select_feature.R --proot ${params.root} \
    --obj ${obj} \
    --iter 5 \
    --resolution 0.2
  """
}

workflow {

  main:
  depth_ch = GetDepth(files(params.frags))
  cn_ch = Annotate(file(params.initobj), depth_ch, files(params.frags))
  raw_cmat_ch = GetRawCount(files(params.bam), file(params.gtf), files(params.cbc))
  meta_ch = GetSplitMeta(
    channel.fromList(['stf_2', 'stf_3', 'stf_4', 'stf_5']).combine(channel.fromList(['class,spatial_origin'])).combine(cn_ch.metadata)
  ).map { m ->
    return [m.baseName, m]
  }
  abam_ch = channel.fromList(['stf_2', 'stf_3', 'stf_4', 'stf_5'])
    .map { l ->
      def bam_path = file("${params.bam}/${l}/outs/atac_possorted_bam.bam")
      return [l, bam_path]
    }
    .combine(meta_ch, by: 0)

  gbam_ch = channel.fromList(['stf_2', 'stf_3', 'stf_4', 'stf_5'])
    .map { l ->
      def bam_path = file("${params.bam}/${l}/outs/gex_possorted_bam.bam")
      return [l, bam_path]
    }
    .combine(meta_ch, by: 0)

  cabam_ch = SplitATACBam(abam_ch)
  cgbam_ch = SplitGexBam(gbam_ch)

  abw_ch = GetATACbw(cabam_ch.flatten())
  gbw_ch = GetGEXbw(cgbam_ch.flatten())

  agbw_ch = AvgGEXbw(
    gbw_ch.map { f ->
      return tuple(f.baseName, f)
    }.groupTuple()
  )
  aabw_ch = AvgATACbw(
    abw_ch.flatten().map { f ->
      def cluster = f.baseName[0..-3]
      return tuple(cluster, f)
    }.groupTuple()
  )

  publish:
  wobj = cn_ch.wobj
  onobj = cn_ch.onobj
  offobj = cn_ch.offobj
  fig_f2 = cn_ch.figs.flatten()
  sup_s2 = cn_ch.supfigs.flatten()
  raw_count = raw_cmat_ch
  atac_bw = aabw_ch
  gex_bw = agbw_ch
}

output {
  wobj {
    path { f ->
      f >> "int/objects/${f.name}"
    }
  }
  onobj {
    path { f ->
      f >> "int/objects/${f.name}"
    }
  }
  offobj {
    path { f ->
      f >> "int/objects/${f.name}"
    }
  }
  fig_f2 {
    path { f ->
      f >> "result/fig/f2/${f.name}"
    }
  }
  sup_s2 {
    path { f ->
      f >> "result/sup_fig/s2/${f.name}"
    }
  }
  raw_count {
    path 'int/'
  }
  atac_bw {
    path 'result/bw/ATAC'
  }
  gex_bw {
    path 'result/bw/GEX'
  }
}
