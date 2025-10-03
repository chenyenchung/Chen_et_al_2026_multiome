params.fastq = '/scratch/cgsb/desplan/Libraries_raw_data/2025-02-27-Azenta-30-1154387723/00_fastq/*.fastq.gz'
params.fasta = '/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
params.fastqp = '/scratch/cgsb/desplan/Libraries_raw_data/2025-02-27-Azenta-30-1154387723/00_fastq/*_R{1,2}_001.fastq.gz'

process fastqc {
  cpus '4'
  module 'fastqc/0.11.9'
  memory '1GB'
  time '1h'
  publishDir 'int/wgs_fastqc/', mode: 'copy', overwrite: false

  input:
  path fastq

  output:
  path "${fastq.baseName}.html"

  script:
  """
  fastqc -t ${task.cpus} ${fastq}
  mv *.html ${fastq.baseName}.html
  """
}

process bwa2_index {
  cpus '1'
  memory '4GB'
  time '10m'

  input:
  path fasta

  output:
  path "dm6.*"

  script:
  """
  bwa-mem2 index -p dm6 ${fasta}
  """
}

process cutadapt {
  cpus '8'
  memory '1GB'
  time '30m'

  module 'cutadapt/4.9'

  input:
  tuple val(sampleId), file(reads)

  output:
  tuple val("${sampleId}"), path('r1.fastq.gz'), path('r2.fastq.gz')

  script:
  """
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
    -o r1.fastq.gz \
    -p r2.fastq.gz \
    -q 20 \
    -j ${task.cpus} \
    ${reads}
  """
}

process bwa2_align {
  cpus '8'
  memory '4GB'
  time '4h'

  input:
  tuple val(sampleId), path(r1), path(r2)
  path idx

  output:
  path "${sampleId}.sam"

  script:
  """
  bwa-mem2 mem -t ${task.cpus} 'dm6' ${r1} ${r2} > ${sampleId}.sam
  """
}

process to_bam {
  cpus '8'
  module 'samtools/intel/1.20'
  memory '16GB'
  time '2h'
  publishDir 'int/wgs_bam/', mode: 'copy', overwrite: false

  input:
  path sam

  output:
  path "${sam.baseName}.bam"
  path "${sam.baseName}.bam.bai"

  script:
  """
  samtools view -b -q 30 -@ ${task.cpus} ${sam} |\
  samtools sort -@ ${task.cpus} -o ${sam.baseName}.bam
  samtools index -@ ${task.cpus} ${sam.baseName}.bam
  """
}

workflow {
  idx = channel.fromPath(params.fasta)
  idxf = bwa2_index(idx).collect()

  fqs = channel.fromPath(params.fastq)
  fastqc(fqs)

  fqp = channel.fromFilePairs(params.fastqp, checkIfExists: true)
  fqt = cutadapt(fqp)
  bwa2_align(fqt, idxf)
    | to_bam
}
