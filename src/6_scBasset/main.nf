#!/usr/bin/env nextflow
nextflow.enable.types = true
nextflow.enable.strict = true

params {
  root: String = '/scratch/ycc520/thesis'
  obj: Path = '/scratch/ycc520/thesis/int/objects/annotated.rds'
  varp: Path = '/scratch/ycc520/thesis/int/hvg_np/varp.rds'
  fasta: Path = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
  sif: Path = '/scratch/ycc520/thesis/static/cuda_11.2.2-cudnn8-devel-ubi8.sif'
  ext3: Path = '/scratch/ycc520/thesis/static/scbasset.ext3'
  prep_src: Path = '/scratch/ycc520/thesis/src/6_scBasset/bin/scbasset_preprocess.py'
  train_src: Path = '/scratch/ycc520/thesis/src/6_scBasset/bin/scbasset_train.py'
}

process ExtractMatrix {
  module 'r/4.5.1'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path
  varp: Path

  output:
  barcodes: Path = file('barcodes.txt')
  peaks: Path = file('peaks.bed')
  matrix: Path = file('peaks.mm')
  meta: Path = file('metadata.csv')

  script:
  """
  extract_mat.R --proot ${params.root} \
    --obj ${obj} \
    --p ${varp}
  """
}

process MakeH5AD {
  conda 'pandas scipy scanpy'
  cpus '1'
  memory '8GB'
  time '30m'

  input:
  _barcodes: Path
  _peaks: Path
  _matrix: Path
  _meta: Path

  output:
  h5ad: Path = file('atac.h5ad')

  script:
  """
  make_h5ad.py
  """
}

process PrepTrain {
  cpus '1'
  memory '8GB'
  time '1h'

  input:
  h5ad: Path
  fasta: Path
  sif: Path
  overlay: Path
  script: Path

  output:
  processed: Path = file('processed', type: 'dir')

  script:
  """
  singularity exec --overlay ${overlay}:ro \
    ${sif} \
    bash -c "source /ext3/bin/enable_mamba.sh; \
             micromamba run -n scBasset python ${script} \
               --ad_file ${h5ad} \
               --input_fasta ${fasta}"
  """
}

process Train {
  cpus '1'
  memory '16GB'
  time '12h'
  clusterOptions '--gres=gpu:1'

  input:
  data: Path
  sif: Path
  overlay: Path
  script: Path

  output:
  model: Set<Path> = files('output')

  script:
  """
  singularity exec --nv --overlay ${overlay}:ro \
    ${sif} \
    bash -c "source /ext3/bin/enable_mamba.sh; \
             micromamba run -n scBasset python ${script} \
               --input_folder ${data} --epochs 1000"
  """
}

workflow {

  main:
  mat_ch = ExtractMatrix(params.obj, params.varp)
  h5ad_ch = MakeH5AD(mat_ch.barcodes, mat_ch.peaks, mat_ch.matrix, mat_ch.meta)
  prep_ch = PrepTrain(
    h5ad_ch,
    params.fasta,
    params.sif,
    params.ext3,
    params.prep_src,
  )
  train_ch = Train(prep_ch, params.sif, params.ext3, params.train_src)

  publish:
  model = train_ch.flatten()
  obj = prep_ch
}

output {
  model {
    path { f ->
      f >> "scBasset/${f.name}"
    }
  }
  obj {
    path "scBasset/"
  }
}
