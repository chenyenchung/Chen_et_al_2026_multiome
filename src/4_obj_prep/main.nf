#!/usr/bin/env nextflow
params.root = '/scratch/ycc520/thesis'
params.bam = '/scratch/ycc520/thesis/data/bam_new'
params.libs = ['stf_2', 'stf_3', 'stf_4', 'stf_5']
params.gtf = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf'
params.fai = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa.fai'
params.demux = '/scratch/ycc520/thesis/int/demux/'
params.freemux = '/scratch/ycc520/thesis/int/freemux/'
params.cbcdir = '/scratch/ycc520/thesis/int/permissive_cbc'
params.mkrds = '/scratch/ycc520/thesis/data/NN_asset/MarkersNNP15.rds'
params.img = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/image'
params.anno = '/scratch/ycc520/thesis/static/ozel_2020_lut.csv'
params.multivi = '/scratch/ycc520/thesis/src/4_obj_prep/bin/run_multivi.py'
params.scvi = '/scratch/ycc520/thesis/src/4_obj_prep/bin/run_scvi.py'
params.peakvi = '/scratch/ycc520/thesis/src/4_obj_prep/bin/run_peakvi.py'
params.blacklist = '/scratch/ycc520/thesis/static/dm6-blacklist.v2.bed.gz'

process GetPrelimPeaks {
  module 'r/4.5.1'
  cpus '1'
  memory '6GB'
  time '30m'

  input:
  path cellranger

  output:
  path "consensus.bed"

  script:
  """
  get_consensus_peaks.R \
    --proot ${params.root}
  """
}

process FilterFragments {
  conda 'bedtools samtools'
  cpus '1'
  memory '6GB'
  time '30m'

  input:
  tuple val(lib), path(peak), path(frag), path(cbc)

  output:
  path "${lib}.tsv.gz"

  script:
  """
  zcat ${frag} |\
    awk -v bc=${cbc} 'BEGIN{OFS="\t"; while((getline b < bc)>0) v[b]=1} \
      !/^#/ && (\$4 in v) {print \$1,\$2,\$2+1,\$4,\$5; print \$1,\$3,\$3+1,\$4,\$5}' |\
    bedtools intersect -a - -b ${peak} -wa -wb |\
    bgzip > ${lib}.tsv.gz

  """
}

process MakeObject {
  module 'r/4.5.1'
  cpus '2'
  memory '32GB'
  time '1h'

  publishDir 'int/objects', mode: 'copy'

  input:
  path cbc
  path gtf
  path data
  path demux
  path freemux
  path intcount

  output:
  path "base_obj.rds"

  script:
  """
  prep_object.R --proot ${params.root} \
    --cbc ${cbc} \
    --gtf ${gtf} \
    --data ${data} \
    --demux ${demux} \
    --freemux ${freemux} \
    --blacklist ${params.blacklist} \
    --out base_obj.rds
  """
}

process FeatureSelection {
  module 'r/4.5.1'
  cpus '1'
  memory '32GB'
  time '3h'

  publishDir 'int/hvg_init', mode: 'copy'

  input:
  path obj

  output:
  tuple path('varg.rds'), path('varp.rds')

  script:
  """
  select_feature.R --proot ${params.root} \
    --obj ${obj} \
    --iter 5 \
    --resolution 0.2
  """
}

process ExtractMat {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'
  publishDir 'int/nn_pred/', mode: 'copy'

  input:
  path obj

  output:
  path "gex_mat.rds"

  script:
  """
  extract_gexmat.R --proot ${params.root} \
    --obj ${obj} \
    --mkrds ${params.mkrds}
  """
}

process NnPredict {
  cpus '2'
  memory '48GB'
  time '30m'
  publishDir 'int/', mode: 'copy'

  input:
  path gex

  output:
  path "nn_pred/"

  script:
  """
  singularity run ${params.img}/nn.sif \
    --stage P15 \
    --matrix ${gex} \
    --output nn_pred/
  """
}

process AnnotateObject {
  module 'r/4.5.1'
  cpus '2'
  memory '32GB'
  time '30m'
  publishDir 'int/objects', mode: 'copy'

  input:
  path obj
  path nnout
  path anno
  tuple path(varg), path(varp)

  output:
  path "base_ann_obj.rds"

  script:
  """
  annotate.R --proot ${params.root} \
    --obj ${obj} \
    --p ${varp} \
    --g ${varg} \
    --resolution 1 \
    --nnout ${nnout} \
    --anno ${anno}
  """
}

process FormatFragment {
  conda 'samtools'
  cpus '1'
  memory '4GB'
  time '30m'

  input:
  tuple val(lib), path(frag)

  output:
  tuple file("${lib}_all.tsv.gz"), file("${lib}_all.tsv.gz.tbi")

  script:
  """
  zcat ${frag} |\
    awk '
  BEGIN{OFS="\t"}
  !/^#/ {print \$1, \$2, \$3, "${lib}#" \$4, \$5}
  ' |\
  bgzip > "${lib}_all.tsv.gz"
  tabix -p bed ${lib}_all.tsv.gz
  """
}

process ManualSubset {
  module 'r/4.5.1'
  cpus '1'
  memory '32GB'
  time '1h'
  publishDir 'int/objects', mode: 'copy'

  input:
  path obj
  path frags

  output:
  path 'intermediate.rds', emit: obj
  path 'metadata.csv', emit: meta

  script:
  """
  manual_subset.R \
    --proot ${params.root} \
    --obj ${obj}
  """
}

process ExtractBarcode {
  module 'r/4.5.1'
  cpus '1'
  memory '24GB'
  time '15m'

  input:
  path obj

  output:
  path 'stf_*_cluster_*.txt'

  script:
  """
  extract_bc.R \
    --proot ${params.root} \
    --obj ${obj}
  """
}

process FilterClusterFragments {
  conda 'samtools'
  cpus '1'
  memory '6GB'
  time '30m'

  input:
  tuple val(lib), val(cluster), file(cbc), file(frag)

  output:
  tuple val("${cluster}"), file("${cbc.baseName}_insertions.tsv.gz")

  script:
  """
  zcat ${frag} |\
    awk -v bc=${cbc} 'BEGIN{OFS="\t"; while((getline b < bc)>0) v[b]=1} \
      !/^#/ && (\$4 in v) {print \$1,\$2,\$2+1,\$4,\$5; print \$1,\$3,\$3+1,\$4,\$5}' |\
    bgzip > ${cbc.baseName}_insertions.tsv.gz

  """
}

process MergeClusterFragments {
  conda 'samtools'
  cpus '1'
  memory '6GB'
  time '30m'

  input:
  tuple val(cluster), file(frags)

  output:
  tuple val("${cluster}"), file("${cluster}_insertions.bed.gz")

  script:
  """
  zcat ${frags} |\
    sort -k1,1 -k2,2n |\
    bgzip > ${cluster}_insertions.bed.gz
  """
}

process ClusterPeak {
  conda 'python=2.7 macs2'
  cpus '1'
  memory '8GB'
  time '1h'

  input:
  tuple val(cluster), file(ins)

  output:
  file "${cluster}_summits.bed"

  script:
  """
  macs2 callpeak \
    -t ${ins} \
    -f BED \
    -g 142573017 \
    --nomodel \
    --shift -75 \
    --extsize 150 \
    -n ${cluster} \
    --keep-dup all \
    --call-summits \
    -q 0.05
  """
}

process GreedyConsensusPeak {
  cpus '1'
  memory '8GB'
  time '30m'
  publishDir 'int/', mode: 'copy'

  input:
  file summits

  output:
  file 'consensus_peaks.bed'

  script:
  """
  cat ${summits} | sort -k5,5rn | awk -v width=150 -v fai=${params.fai} '
  BEGIN {
    OFS="\t"
    while ((getline < fai) > 0) chrlen[\$1] = \$2
  }
  {
      chr = \$1; pos = \$2

      # Skip if would extend past chromosome ends
      if (pos - width < 0 || pos + width > chrlen[chr]) next

      bin = int(pos / width)
      dominated = 0
      
      for (b = bin - 1; b <= bin + 1; b++) {
          key = chr":"b
          if (key in kept) {
              d = pos - kept[key]
              if (d < 0) d = -d
              if (d < width) { dominated = 1; break }
          }
      }
      
      if (!dominated) {
          kept[chr":"bin] = pos
          print chr, pos - width, pos + width, \$4, \$5
      }
  }' | sort -k1,1 -k2,2n > consensus_peaks.bed
  """
}

process PileBarcode {
  conda 'bedtools samtools'
  cpus '1'
  memory '6GB'
  time '30m'
  publishDir 'int/filtered_barcode', mode: 'copy'

  input:
  tuple val(lib), file(cbcs)

  output:
  tuple val("${lib}"), file("${lib}.txt")

  script:
  """
  cat ${cbcs} | sort > "${lib}.txt"
  """
}

process QuantifyFragments {
  conda 'bedtools samtools'
  cpus '1'
  memory '6GB'
  time '30m'

  input:
  tuple val(lib), path(cbc), path(peak), path(frag)

  output:
  path "${lib}.tsv.gz"

  script:
  """
  zcat ${frag} |\
    awk -v bc=${cbc} 'BEGIN{OFS="\t"; while((getline b < bc)>0) v[b]=1} \
      !/^#/ && (\$4 in v) {print \$1,\$2,\$2+1,\$4,\$5; print \$1,\$3,\$3+1,\$4,\$5}' |\
    bedtools intersect -a - -b ${peak} -wa -wb |\
    bgzip > ${lib}.tsv.gz
  """
}

process MakeIntermediateObject {
  module 'r/4.5.1'
  cpus '1'
  memory '32GB'
  time '1h'

  publishDir 'int/objects', mode: 'copy'

  input:
  path cbc
  path gtf
  path data
  path demux
  path freemux
  path intcount
  path meta

  output:
  path "filtered_obj.rds"

  script:
  """
  prep_object.R --proot ${params.root} \
    --cbc ./ \
    --gtf ${gtf} \
    --data ${data} \
    --demux ${demux} \
    --freemux ${freemux} \
    --blacklist ${params.blacklist} \
    --meta ${meta} \
    --out filtered_obj.rds
  """
}

process FeatureSelectionNP {
  module 'r/4.5.1'
  cpus '1'
  memory '36GB'
  time '4h'

  publishDir 'int/hvg_np', mode: 'copy'

  input:
  path obj

  output:
  tuple path('varg.rds'), path('varp.rds')

  script:
  """
  select_feature.R --proot ${params.root} \
    --obj ${obj} \
    --iter 5 \
    --resolution 0.2
  """
}

process MakeInsertionFileNP {
  conda 'samtools'
  cpus '1'
  memory '6GB'
  time '30m'
  publishDir 'int/insertion_np', mode: 'copy'

  input:
  tuple val(lib), path(frag), path(cbc)

  output:
  path "${lib}.tsv.gz"
  path "${lib}.tsv.gz.tbi"

  script:
  """
  zcat ${frag} |\
    awk -v bc=${cbc} 'BEGIN{OFS="\t"; while((getline b < bc)>0) v[b]=1} \
      !/^#/ && (\$4 in v) {print \$1,\$2,\$2+1,"${lib}#" \$4,\$5; print \$1,\$3,\$3+1,"${lib}#" \$4,_\$5}' |\
      sort -k1,1 -k2,2n |\
    bgzip > ${lib}.tsv.gz
  tabix -p bed ${lib}.tsv.gz
  """
}

process MakeWNNObject {
  module 'r/4.5.1'
  cpus '2'
  memory '36GB'
  time '1h'
  publishDir 'int/objects', mode: 'copy'

  input:
  path obj
  tuple path(varg), path(varp)

  output:
  path 'filtered_wnn.rds'

  script:
  """
  make_wnn_obj.R --proot ${params.root} \
    --obj ${obj} \
    --p ${varp} \
    --g ${varg} \
    --resolution 0.5 
  """
}

process MatrixExtSCVI {
  module 'r/4.5.1'
  cpus '1'
  memory '12GB'
  time '30m'

  input:
  path obj
  tuple path(varg), path(varp)

  output:
  tuple file("genes.txt"), file("barcodes.txt"), file("gex.mm"), file("peaks.txt"), file("peaks.mm"), file("metadata.csv")

  script:
  """
  mat_ext_scvi.R --proot ${params.root} \
    --obj ${obj} \
    --p ${varp} \
    --g ${varg}
  """
}

process MakeMudata {
  cpus '1'
  memory '8GB'
  time '30m'
  conda 'muon'

  input:
  tuple file(genes), file(barcodes), file(gex_mat), file(peaks), file(peak_mat), file(meatdata)

  output:
  file 'obj.h5mu'

  script:
  """
  make_mudata.py
  """
}

process RunMultiVI {
  cpus '1'
  memory '16GB'
  time '2h'
  clusterOptions '--gres=gpu:1'
  publishDir 'int/scvi_tools/', mode: 'copy'

  input:
  tuple path(mudata), path(runscript), val(k)

  output:
  path "multiVI/*"

  script:
  """
  DEP_PREFIX="/projects/rps/cgsb/desplan/File_exchange/Yen_ref"
  IMAGE_PREFIX="\${DEP_PREFIX}/image"
  OVERLAY_PREFIX="\${DEP_PREFIX}/overlay"

  singularity exec --nv --overlay "\${OVERLAY_PREFIX}/scvi_nv.ext3:ro" \
    "\${IMAGE_PREFIX}/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif" \
    bash -c "source /ext3/bin/enable_micromamba.sh; \
      micromamba run -n flyvis python ${runscript} \
        --h5mu ${mudata} \
        --nfactor ${k} \
        --batch library \
        --penalty Jeffreys"
  """
}

process RunSCVI {
  cpus '1'
  memory '16GB'
  time '1h'
  clusterOptions '--gres=gpu:1'
  publishDir 'int/scvi_tools/', mode: 'copy'

  input:
  tuple path(mudata), path(runscript), val(k)

  output:
  path "scVI/*"

  script:
  """
  DEP_PREFIX="/projects/rps/cgsb/desplan/File_exchange/Yen_ref"
  IMAGE_PREFIX="\${DEP_PREFIX}/image"
  OVERLAY_PREFIX="\${DEP_PREFIX}/overlay"

  singularity exec --nv --overlay "\${OVERLAY_PREFIX}/scvi_nv.ext3:ro" \
    "\${IMAGE_PREFIX}/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif" \
    bash -c "source /ext3/bin/enable_micromamba.sh; \
      micromamba run -n flyvis python ${runscript} \
        --h5mu ${mudata} \
        --nfactor ${k} \
        --batch library"
  """
}

process RunPeakVI {
  cpus '1'
  memory '16GB'
  time '1h'
  clusterOptions '--gres=gpu:1'
  publishDir 'int/scvi_tools/', mode: 'copy'

  input:
  tuple path(mudata), path(runscript), val(k)

  output:
  path "peakVI/*"

  script:
  """
  DEP_PREFIX="/projects/rps/cgsb/desplan/File_exchange/Yen_ref"
  IMAGE_PREFIX="\${DEP_PREFIX}/image"
  OVERLAY_PREFIX="\${DEP_PREFIX}/overlay"

  singularity exec --nv --overlay "\${OVERLAY_PREFIX}/scvi_nv.ext3:ro" \
    "\${IMAGE_PREFIX}/cuda12.8.1-cudnn9.8.0-ubuntu24.04.2.sif" \
    bash -c "source /ext3/bin/enable_micromamba.sh; \
      micromamba run -n flyvis python ${runscript} \
        --h5mu ${mudata} \
        --nfactor ${k} \
        --batch library"
  """
}

process SelectDim {
  module 'r/4.5.1'
  cpus '16'
  memory '12GB'
  time '2h'
  publishDir 'result/sup_fig/X1_select_dimensions', mode: 'copy'

  input:
  path obj

  output:
  path 'gex_dim_jaccard.csv'
  path 'PCA_jaccard_heatmap.png'
  path 'PCA_jaccard_curve.png'
  path 'atac_dim_jaccard.csv'
  path 'LSI_jaccard_heatmap.png'
  path 'LSI_jaccard_curve.png'

  script:
  """
  select_dim.R --proot ${params.root} \
    --obj ${obj} \
    --ncpus ${task.cpus}
  """
}

process Cluster {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '1h'
  publishDir './', mode: 'copy', saveAs: { fn ->
    fn.endsWith(".pdf") ? "result/sup_fig/X2_clustering/${fn}" : "int/objects/${fn}"
  }

  input:
  path obj

  output:
  path 'raw_cluster.rds'
  path 'graph_lapacian_eig.pdf'
  path 'graph_lapacian_eig_gap.pdf'
  path 'graph_lapacian_silhouette.pdf'

  script:
  """
  clustering.R --proot ${params.root} \
    --obj ${obj} \
    --resolution 2
  """
}

workflow {
  peak_ch = GetPrelimPeaks(params.bam)
  lib_ch = channel.fromList(params.libs).combine(peak_ch)
    | map { lib, peak ->
      def frag = file("data/bam_new/${lib}/outs/atac_fragments.tsv.gz")
      def cbc = file("int/permissive_cbc/${lib}.txt")
      return [
        lib,
        peak,
        frag,
        cbc,
      ]
    }

  premat_ch = FilterFragments(lib_ch).collect()

  obj_ch = MakeObject(
    file(params.cbcdir),
    file(params.gtf),
    file(params.bam),
    file(params.demux),
    file(params.freemux),
    premat_ch,
  )

  feat_ch = FeatureSelection(obj_ch)
  nn_ch = ExtractMat(obj_ch) | NnPredict
  annobj_ch = AnnotateObject(obj_ch, nn_ch, file(params.anno), feat_ch)

  raw_frag_ch = FormatFragment(
    channel.fromList(params.libs) | map { lib ->
      def frag = file("data/bam_new/${lib}/outs/atac_fragments.tsv.gz")
      return [
        lib,
        frag,
      ]
    }
  ).collect()

  subset_ch = ManualSubset(annobj_ch, raw_frag_ch)
  clusterbc_ch = ExtractBarcode(subset_ch.obj)
    .flatten()
    .map { e ->
      def fname = e.getName()
      def m = (fname =~ /^(.+?)_(cluster_\d+)\.txt$/)
      def frag = "data/bam_new/${m[0][1]}/outs/atac_fragments.tsv.gz"
      return [m[0][1], m[0][2], file(e), file(frag)]
    }

  insertion_ch = FilterClusterFragments(clusterbc_ch)
    | groupTuple(by: 0)
    | MergeClusterFragments

  peak_ch = ClusterPeak(insertion_ch)
    | collect
    | GreedyConsensusPeak

  pbc_ch = clusterbc_ch
    | map { lib, _clust, cbc_file, _frag ->
      return [lib, cbc_file]
    }
    | groupTuple(by: 0)
    | PileBarcode

  pcount_ch = pbc_ch
    | combine(peak_ch)
    | map { lib, clust, peak ->
      def frag = file("data/bam_new/${lib}/outs/atac_fragments.tsv.gz")
      return [lib, clust, peak, frag]
    }
    | QuantifyFragments

  nobj_ch = MakeIntermediateObject(
    pbc_ch.map { _lib, cbc ->
      return cbc
    }.collect(),
    file(params.gtf),
    file(params.bam),
    file(params.demux),
    file(params.freemux),
    pcount_ch.collect(),
    subset_ch.meta,
  )

  newvar_ch = FeatureSelectionNP(nobj_ch)
  npins_ch = MakeInsertionFileNP(
    channel.fromList(params.libs) | map { lib ->
      def frag = file("data/bam_new/${lib}/outs/atac_fragments.tsv.gz")
      return [
        lib,
        frag,
      ]
    } | join(pbc_ch)
  )
  wnnobj_ch = MakeWNNObject(nobj_ch, newvar_ch)
  sel_dim_ch = SelectDim(wnnobj_ch)
  preobj_ch = Cluster(wnnobj_ch)
}
