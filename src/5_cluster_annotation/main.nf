#!/usr/bin/env nextflow
nextflow.preview.types = true
nextflow.enable.strict = true

params.root = '/scratch/ycc520/thesis'
params.initobj = '/scratch/ycc520/thesis/int/objects/raw_cluster.rds'
params.frags = '/scratch/ycc520/thesis/int/insertion_np'

process Annotate {
  module 'r/gcc/4.5.0'
  cpus '1'
  memory '16GB'
  time '30m'

  input:
  obj: Path
  _frags: List<Path>

  output:
  wobj: Path = file('obj/annotated.rds')
  onobj: Path = file('obj/non_obj.rds')
  offobj: Path = file('obj/noff_obj.rds')
  figs: Set<Path> = files('fig/**')
  supfigs: Set<Path> = files('sup/**')

  script:
  """
  annotate.R --proot ${params.root} \
    --obj ${obj}
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
  cn_ch = Annotate(file(params.initobj), files(params.frags))

  publish:
  wobj = cn_ch.wobj
  onobj = cn_ch.onobj
  offobj = cn_ch.offobj
  fig_f2 = cn_ch.figs.flatten()
  sup_s2 = cn_ch.supfigs.flatten()
}

output {
  wobj {
    path { f ->
      f >> "int/objects/${f.name}"
    }
  }
  onobj {
    path { f ->
      f >> "int/${f.name}"
    }
  }
  offobj {
    path { f ->
      f >> "int/${f.name}"
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
}
