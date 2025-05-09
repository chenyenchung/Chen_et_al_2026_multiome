#!/usr/bin/env nextflow
nextflow.preview.output = true

params.lines = '/scratch/ycc520/thesis/data/line_used.txt'
params.input = '/vast/ycc520/data/dgrp/dgrp2_dm6.vcf'
params.gtf = '/vast/ycc520/data/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.88.gtf'
params.driver = '/scratch/ycc520/thesis/int/wgs_variant/gatk_snp.recalibrated.vcf.gz'
params.extern = '/scratch/ycc520/thesis/extern'
params.libpaths = '/scratch/ycc520/thesis/data/bam_new'
params.popscle = '/vast/ycc520/data/popscle.sif'
params.libs = ['stf_2', 'stf_3', 'stf_4', 'stf_5']

process GetDGRPInfo {
  module "bcftools/intel/1.14:samtools/intel/1.14:vcftools/intel/0.1.16"
  cpus '1'
  memory '1GB'
  time '15m'

  input:
  path vcf
  val lines

  output:
  path "lines_used.vcf.gz"

  script:
  def linesf = lines.join(',')
  def max_lines = 1
  """
  source "${params.extern}/popscle_helper_tools/filter_vcf_file_for_popscle.sh"
  subset_samples_from_vcf ${linesf} ${vcf} |\
  only_keep_snps |\
  filter_out_mutations_missing_genotype_for_one_or_more_samples |\
  filter_out_mutations_homozygous_reference_in_all_samples |\
  filter_out_mutations_homozygous_in_all_samples |\
  filter_out_mutations_heterozygous_for_one_or_more_samples |\
  only_keep_mutations_homozygous_in_k_samples ${max_lines} | \
  bgzip > "lines_used.vcf.gz"
  bcftools index "lines_used.vcf.gz"
  """
}

process BcftoolsIntersect {
  cpus '1'
  memory '1GB'
  time '15m'
  module 'bcftools/intel/1.14'

  input:
  path vcf1
  path vcf2

  output:
  tuple path('driver_filtered.vcf.gz'), path('dgrp_filtered.vcf.gz'), emit: shared
  path 'driver_only.vcf.gz', emit: driver_only
  path 'dgrp_asis.vcf.gz', emit: dgrp_asis

  script:
  """
  bcftools index ${vcf1}
  bcftools index ${vcf2}
  bcftools isec -p output \
  -Oz ${vcf1} ${vcf2}
  mv output/0000.vcf.gz driver_only.vcf.gz
  mv output/0001.vcf.gz dgrp_asis.vcf.gz
  mv output/0002.vcf.gz driver_filtered.vcf.gz
  mv output/0003.vcf.gz dgrp_filtered.vcf.gz
  """
}

process GenerateSharedVCF {
  module 'bcftools/intel/1.14:r/gcc/4.4.0'

  cpus '1'
  time '15m'
  memory '1GB'

  input:
  val lib
  path meta
  tuple path(driver), path(dgrp)

  output:
  tuple val("${lib}"), path("${lib}_shared.vcf.gz"), path("${lib}_shared.vcf.gz.csi")

  script:
  """
    generate_synth.R \
      --lib ${lib} \
      --meta ${meta} \
      --dgrp ${dgrp} \
      --driver ${driver} \
      --out "${lib}_shared.vcf.gz"
    gunzip "${lib}_shared.vcf.gz"
    bgzip "${lib}_shared.vcf"
    bcftools index "${lib}_shared.vcf.gz"
    """
}

process GenerateASISVCF {
  module 'bcftools/intel/1.14:r/gcc/4.4.0'
  cpus '1'
  time '15m'
  memory '2GB'

  input:
  val lib
  path meta
  path dgrp

  output:
  tuple val("${lib}"), path("${lib}_asis.vcf.gz"), path("${lib}_asis.vcf.gz.csi")

  script:
  """
  reformat.R \
    --lib ${lib} \
    --meta ${meta} \
    --dgrp ${dgrp} \
    --out "${lib}_asis.vcf.gz"
  gunzip "${lib}_asis.vcf.gz"
  bgzip "${lib}_asis.vcf"
  bcftools index "${lib}_asis.vcf.gz"
  """
}

process MaskHighlyExpGenes {
  module 'r/gcc/4.4.0'
  cpus '1'
  memory '8GB'
  time '15m'

  input:
  val lib
  path h5path
  path gtf

  output:
  tuple val("${lib}"), path("${lib}_mask.bed")

  script:
  """
    get_top.R \
      --lib ${lib} \
      --h5 ${h5path} \
      --gtf ${gtf}
    """
}

process ConcatVCF {
  module "bcftools/intel/1.14:samtools/intel/1.14:vcftools/intel/0.1.16"
  cpus '1'
  memory '1GB'
  time '15m'

  input:
  tuple val(lib), path(shared), path(sharedi), path(dgrp), path(dgrpi), path(mask)

  output:
  tuple val("${lib}"), path("${lib}_info.vcf.gz")

  script:
  """
  source "${params.extern}/popscle_helper_tools/filter_vcf_file_for_popscle.sh"
  bcftools concat -a \
    ${shared} ${dgrp} |\
    bcftools view -T ^${mask} |\
    calculate_AF_AC_AN_values_based_on_genotype_info |\
    bgzip > "${lib}_info.vcf.gz"
    bcftools index "${lib}_info.vcf.gz"
  """
}

process PermissiveFilterBarcodes {
  module 'r/gcc/4.4.0'
  cpus '1'
  time '15m'
  memory '8GB'

  input:
  val lib
  path h5

  output:
  tuple val("${lib}"), path("${lib}.txt")

  script:
  """
  permissive_filter.R \
    --lib ${lib} \
    --h5path ${h5}
  """
}

process PreFilterBAM {
  module 'samtools/intel/1.14:bcftools/intel/1.14:bedtools/intel/2.29.2'
  cpus '1'
  memory '500MB'
  time '3h'

  input:
  tuple val(lib), path(bam_gex), path(bam_atac), path(bclist), path(vcf)

  output:
  tuple val("${lib}"), path("${lib}_gex.bam"), path("${lib}_atac.bam")

  script:
  """
  ${params.extern}/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    ${bam_gex} \
    ${bclist} \
    ${vcf} \
    "${lib}_gex.bam"
  
  ${params.extern}/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh \
    ${bam_atac} \
    ${bclist} \
    ${vcf} \
    "${lib}_atac.bam"
  """
}

process MergeBAM {
  module 'samtools/intel/1.14'
  cpus '8'
  memory '2G'
  time '1h'

  input:
  tuple val(lib), path(gex_bam), path(atac_bam)

  output:
  tuple val("${lib}"), path("${lib}.bam")

  script:
  """
  samtools merge \
    -@ ${task.cpus} \
    -o "${lib}.bam" \
    ${gex_bam} ${atac_bam}
  """
}

process SortVCF {
  module 'bcftools/intel/1.14:samtools/intel/1.14'
  cpus '1'
  memory '2G'
  time '15m'

  input:
  tuple val(lib), path(vcf), path(bam)

  output:
  tuple val("${lib}"), path("${lib}_sorted.vcf.gz")

  script:
  """
  ${params.extern}/popscle_helper_tools/sort_vcf_same_as_bam.sh \
    ${bam} \
    ${vcf} > "${lib}_sorted.vcf.gz"
  """
}

process Mpileup {
  module 'samtools/intel/1.14:bcftools/intel/1.14:bedtools/intel/2.29.2:singularity/3.7.4'
  cpus '1'
  memory '4GB'
  time '6h'

  input:
  tuple val(lib), path(vcf), path(bam), path(bclist)
  path image

  output:
  tuple val("${lib}"), path("${bclist.baseName}_demux/*")

  script:
  """
  mkdir -p "${bclist.baseName}_demux/"
  singularity run ${image} "dsc-pileup \
    --sam ${bam} \
    --vcf ${vcf} \
    --group-list ${bclist} \
    --out ${bclist.baseName}_demux/samples_to_demultiplex.pileup"
  """
}

process Demuxlet {
  module 'samtools/intel/1.14:bcftools/intel/1.14:bedtools/intel/2.29.2:singularity/3.7.4'
  cpus '1'
  memory '4GB'
  time '15m'

  input:
  tuple val(lib), path(plp), path(vcf)
  path image

  output:
  tuple val("${lib}"), path("${plp}.best")

  script:
  """
  singularity run ${image} "demuxlet \
    --plp ${plp}/samples_to_demultiplex.pileup \
    --vcf ${vcf} \
    --field GT \
    --alpha 0 --alpha 0.5 \
    --out ${plp}"
  """
}

process CollectTable {

  cpus '1'
  memory '1GB'
  time '15m'

  input:
  tuple val(lib), path(tbl)

  output:
  path "${lib}.best"

  script:
  """
  awk 'NR == 1 || FNR != 1' ${tbl} \
    > "${lib}.best"
  """
}

workflow {
  main:
  lib_ch = Channel.fromList(params.libs)
  // Get a list of lines that were ever used for subsetting
  Channel.fromPath(params.lines)
    | splitCsv(header: false)
    | map { row -> row[1] }
    | unique
    | collect
    | set { lineid_ch }
  dgrp_ch = GetDGRPInfo(file(params.input), lineid_ch)

  // Find DGRP SNPs that are confounded by driver line genotypes
  driver_ch = BcftoolsIntersect(file(params.driver), dgrp_ch)

  // Generate predicted genotype VCFs
  shared_ch = GenerateSharedVCF(lib_ch, file(params.lines), driver_ch.shared)
  asis_ch = GenerateASISVCF(lib_ch, file(params.lines), driver_ch.dgrp_asis)

  // Mask highly expressed genes from informative variants
  mask_ch = MaskHighlyExpGenes(lib_ch, file(params.libpaths), file(params.gtf))

  // Concatenate VCF files with masking
  join_ch = shared_ch.join(asis_ch).join(mask_ch)
  info_ch = ConcatVCF(join_ch)

  // Permissively (> 300 fragments & > 200 genes) filter each library
  bclist_ch = PermissiveFilterBarcodes(lib_ch, file(params.libpaths))
  demux_vch = bclist_ch.join(info_ch)

  lib_ch
    | map { lib ->
      tuple(
        lib,
        file([params.libpaths, lib, 'outs', 'gex_possorted_bam.bam'].join('/')),
        file([params.libpaths, lib, 'outs', 'atac_possorted_bam.bam'].join('/')),
      )
    }
    | join(demux_vch)
    | set { rawbam_ch }

  fbam_ch = PreFilterBAM(rawbam_ch)
  mbam_ch = MergeBAM(fbam_ch)
  sinfo_ch = SortVCF(info_ch.join(mbam_ch))

  // Run mpileup per 1k barcode to run in parallel
  sbclist_ch = bclist_ch
    | map { lib, bcpath ->
      def chunk = bcpath.splitText(by: 500, file: true)
      return [lib, chunk]
    }
    | transpose
  plp_ch = Mpileup(sinfo_ch.join(mbam_ch).combine(sbclist_ch, by: 0), file(params.popscle))

  vplp_ch = plp_ch
    | map { l, p ->
      def pbase = p[0].getParent()
      return [l, pbase]
    }
    | combine(sinfo_ch, by: 0)

  idemux_ch = Demuxlet(vplp_ch, file(params.popscle))
  out_ch = CollectTable(idemux_ch.groupTuple())

  publish:
  tbls = out_ch
}

output {
    tbls {
    path "demuxlet/"
  }
}
