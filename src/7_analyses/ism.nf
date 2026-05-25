#!/usr/bin/env nextflow
nextflow.enable.types = true
nextflow.enable.strict = true

params {
  root: String = '/scratch/ycc520/thesis'
  bed: Path = '/scratch/ycc520/thesis/static/reporter_fragments.bed'
  fasta: Path = '/projects/rps/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
  motif: Path = '/scratch/ycc520/thesis/static/CisBPDrosophilaALL.meme'
  anndata: Path = '/scratch/ycc520/thesis/result/scBasset/processed/ad.h5ad'
  model: Path = '/scratch/ycc520/thesis/result/scBasset/output/best_model.h5'
  sif: Path = '/scratch/ycc520/thesis/static/cuda_11.2.2-cudnn8-devel-ubi8.sif'
  ext3: Path = '/scratch/ycc520/thesis/static/scbasset.ext3'
  ism_src: Path = '/scratch/ycc520/thesis/src/7_analyses/bin/scbasset_ism.py'
  groupby_pairs: String = 'spatial_origin,class;temporal_identity,class'
  seq_len: Integer = 1344
  fimo_thres: Float = 1e-3
  mutation_batch_size: Integer = 256
  top_windows: Integer = 100
  window_size: Integer = 10
}

process PrepareISMRegions {
  cpus '1'
  memory '4GB'
  time '30m'

  input:
  bed: Path
  fasta: Path
  sif: Path
  overlay: Path
  script_file: Path

  output:
  prepared: Path = file('prepared', type: 'dir')

  script:
  """
  singularity exec --overlay ${overlay}:ro \
    ${sif} \
    bash -c "source /ext3/bin/enable_mamba.sh; \
             micromamba run -n scBasset python ${script_file} prepare \
               --bed ${bed} \
               --fasta ${fasta} \
               --outdir prepared \
               --seq-len ${params.seq_len}"
  """
}

process RunFIMO {
  conda 'meme bedtools'
  cpus '1'
  memory '8GB'
  time '1h'

  input:
  prepared: Path
  motif: Path

  output:
  fimo: Path = file('fimo_out', type: 'dir')

  script:
  """
  fimo --verbosity 1 \
    --thresh ${params.fimo_thres} \
    --oc fimo_out \
    ${motif} \
    ${prepared}/regions.fa
  """
}

process RunScBassetISM {
  cpus '1'
  memory '32GB'
  time '1h'
  clusterOptions '--gres=gpu:1'

  input:
  tuple(pair_id: String, pair: String, prepared: Path, fimo: Path)
  sif: Path
  overlay: Path
  script_file: Path
  anndata: Path
  model: Path
  motif: Path

  output:
  outdir: Path = file("${pair_id}", type: 'dir')

  script:
  """
  mkdir -p "${pair_id}"
  singularity exec --nv --overlay ${overlay}:ro \
    ${sif} \
    bash -c "source /ext3/bin/enable_mamba.sh; \
             micromamba run -n scBasset python ${script_file} score \
               --regions-tsv ${prepared}/regions.tsv \
               --fasta ${prepared}/regions.fa \
               --fimo-tsv ${fimo}/fimo.tsv \
               --meme ${motif} \
               --anndata ${anndata} \
               --model ${model} \
               --groupby-pair '${pair}' \
               --outdir ${pair_id} \
               --seq-len ${params.seq_len} \
               --mutation-batch-size ${params.mutation_batch_size} \
               --top-windows ${params.top_windows} \
               --window-size ${params.window_size}"
  """
}

workflow {

  main:
  prepared_ch = PrepareISMRegions(
    params.bed,
    params.fasta,
    params.sif,
    params.ext3,
    params.ism_src,
  )
  fimo_ch = RunFIMO(prepared_ch, params.motif)

  pair_ch = channel.fromList(params.groupby_pairs.split(';').findAll { it.trim() })
    .map { pair ->
      def cleaned = pair.trim()
      def pair_id = cleaned
        .replace(',', '__')
        .replaceAll(/[^A-Za-z0-9_.-]+/, '_')
        .replaceAll(/^_+|_+$/, '')
      return tuple(pair_id, cleaned)
    }

  score_input_ch = pair_ch
    .combine(prepared_ch)
    .combine(fimo_ch)
    .map { item ->
      return tuple(item[0], item[1], item[2], item[3])
    }

  ism_ch = RunScBassetISM(
    score_input_ch,
    params.sif,
    params.ext3,
    params.ism_src,
    params.anndata,
    params.model,
    params.motif,
  )

  publish:
  prepared = prepared_ch
  fimo = fimo_ch
  ism = ism_ch
}

output {
  prepared {
    path "scBasset_ism/prepared"
  }
  fimo {
    path "scBasset_ism/fimo"
  }
  ism {
    path "scBasset_ism"
  }
}
