params.ref='/vast/ycc520/data/ens88/'

process arc_align {
  cpus '24'
  memory '32GB'
  time '18h'
  publishDir 'data/bam_new', mode: 'copy', overwrite: false

  input:
    path sampleSheet
    path ref

  output:
    path "${sampleSheet.baseName}/*"

  script:
  """
  cellranger-arc count --id=${sampleSheet.baseName} \
       --reference=${ref} \
       --libraries=${sampleSheet} \
       --localcores=${task.cpus} \
       --disable-ui
  """
}

workflow {
  sheets = channel.fromPath("data/stf_*.csv")
  arc_align(sheets, file(params.ref))
}
