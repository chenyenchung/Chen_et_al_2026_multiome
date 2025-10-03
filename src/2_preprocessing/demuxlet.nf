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


process GenerateASISVCF {
  module 'bcftools/intel/1.14:r/gcc/4.5.0'
  cpus '1'
  time '15m'
  memory '2GB'

  input:
  val lib
  path meta
  path dgrp

  output:
  tuple val("${lib}"), path("${lib}_asis.vcf.gz")

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
  module 'bcftools/intel/1.14:r/gcc/4.5.0'
  cpus '1'
  memory '8GB'
  time '15m'

  input:
  tuple val(lib), path(vcf)
  path h5path
  path gtf

  output:
  tuple val("${lib}"), path("${lib}_masked.vcf.gz")

  script:
  """
  get_top.R \
    --lib ${lib} \
    --h5 ${h5path} \
    --gtf ${gtf}
  bcftools view ${vcf} \
    -Oz -o "${lib}_masked.vcf.gz" \
    -T "^${lib}_mask.bed"
  """
}

process PermissiveFilterBarcodes {
  publishDir 'int/permissive_cbc/', mode: 'copy'
  module 'r/gcc/4.5.0'
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

process FilterByCoverage {
  module 'bedops/intel/2.4.39:bcftools/intel/1.14:bedtools/intel/2.29.2:samtools/intel/1.14'

  input:
  tuple val(lib), path(bam_gex), path(bam_atac), path(vcf)

  output:
  tuple val("${lib}"), path("${lib}_info.vcf.gz")

  script:
  """
  gunzip -f ${vcf}
  mv *.vcf input.vcf
  
  vcf2bed < input.vcf |\
    awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$4}' |\
    sort -k1,1 -k2,2n > dgrp.bed
  bedtools coverage -a dgrp.bed \
    -b ${bam_gex} -counts -sorted |\
    awk '\$5 >= 10' > coverage.bed

  bedtools coverage -a dgrp.bed \
    -b ${bam_atac} -counts -sorted |\
    awk '\$5 >= 10' >> coverage.bed

  sort -k1,1 -k2,2n coverage.bed |\
    bedtools merge -i - > coverage_merged.bed

  bgzip input.vcf
  bcftools index input.vcf.gz
  source "${params.extern}/popscle_helper_tools/filter_vcf_file_for_popscle.sh"
  bcftools view -Oz -R coverage_merged.bed \
    input.vcf.gz |\
    calculate_AF_AC_AN_values_based_on_genotype_info > ${lib}_info.vcf.gz
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
  cpus '1'
  memory '2GB'
  time '4h'

  input:
  tuple val(lib), path(vcf), path(bam), path(bclist)
  path image

  output:
  tuple val("${lib}"), path("${bclist.baseName}.pileup*")

  script:
  """
  mkdir -p "${bclist.baseName}_demux/"
  singularity run ${image} "dsc-pileup \
    --sam ${bam} \
    --vcf ${vcf} \
    --group-list ${bclist} \
    --out ${bclist.baseName}.pileup"
  """
}

process IntegratePL {
  publishDir 'int/mpileup/', mode: 'copy'
  cpus '1'
  memory '4GB'
  time '2h'
  conda 'python=3.12'

  input:
  tuple val(lib), path(demux_slice)

  output:
  tuple val("${lib}"), path("${lib}.pileup.*.gz")

  script:
  """
#!/usr/bin/env python
import os
import gzip
import re
import hashlib
import shutil

lib = \'${lib}\'
files = os.listdir('./')
prefixes = []
for f in files:
  if re.search(r'\\.gz\$', f):
    p = re.sub(r'\\.pileup\\..{3}\\.gz', '', f)
    if p not in prefixes:
      prefixes.append(p)

uid = 0
uid_dict= {}
cel_header = True
plp_header = True

var_path = []
for d in prefixes:
    ## Deal with .cel files
    with gzip.open(f'./{d}.pileup.cel.gz', 'rt') as f:
        with gzip.open(os.path.join('./', f'{lib}.pileup.cel.gz'), 'at') as o:
            lc = 0
            for line in f:
                if lc == 0:
                    if cel_header:
                        o.write(line)
                        cel_header = False
                else:
                    content = line.strip('\\n')
                    content = content.split('\\t')
                    # Store the barcode IDs (first field) to a dictionary
                    uid_dict[str(content[0])] = str(uid)
                    # and recode incrementally across slices
                    content[0] = str(uid)
                    uid = uid + 1
                    out = '\\t'.join(content) + '\\n'
                    o.write(out)
                lc = lc + 1

    ## Deal with .umi files
    with gzip.open(f'./{d}.pileup.umi.gz', 'rt') as f:
        with gzip.open(os.path.join('./', f'{lib}.pileup.umi.gz'), 'at') as o:
            ulc = 0
            for line in f:
                content = line.strip('\\n').split('\\t')
                # Recode barcode IDs with existing dictionary
                content[0] = uid_dict[content[0]]
                o.write('\\t'.join(content) + '\\n')
                ulc = ulc + 1

    ## Deal with .plp files
    with gzip.open(f'./{d}.pileup.plp.gz', 'rt') as f:
        with gzip.open(os.path.join('./', f'{lib}.pileup.plp.gz'), 'at') as o:
            plc = 0
            for line in f:
                if plc == 0:
                    if plp_header:
                        o.write(line)
                        plp_header = False
                else:
                    content = line.strip('\\n')
                    content = content.split('\\t')
                    # Recode barcode IDs with existing dictionary
                    content[0] = uid_dict[content[0]]
                    out = '\\t'.join(content) + '\\n'
                    o.write(out)
                plc = plc + 1

    var_path.append(f'./{d}.pileup.var.gz')

var_hash = dict()
for v in var_path:
    with open(v, 'rb') as f:
        hash = hashlib.file_digest(f, 'sha256').hexdigest()
        if hash in var_hash:
            var_hash[hash] += 1
        else:
            var_hash[hash] = 1


if len(var_hash) == 1:
    shutil.copyfile(var_path[0], f'{lib}.pileup.var.gz')
else:
    raise ValueError(f'In {lib}, there is more than one version of .var.gz.')

print(f'[INFO] Completed integrating pileup data from {uid} cells.')
"""
}

process Demuxlet {
  cpus '1'
  memory '32GB'
  time '30m'
  publishDir 'int/demux/', mode: 'copy'

  input:
  tuple val(lib), path(plp), path(vcf)
  path image

  output:
  tuple val("${lib}"), path("${lib}.best")

  script:
  """
  singularity run ${image} "demuxlet \
    --plp ${plp}/${lib}.pileup \
    --vcf ${vcf} \
    --field GT \
    --alpha 0 --alpha 0.5 \
    --out ${lib}"
  """
}

process Freemux {
  cpus '1'
  memory '48GB'
  time '8h'
  publishDir 'int/freemux/', mode: 'copy'

  input:
  tuple val(lib), path(plp)
  path image

  output:
  tuple val("${lib}"), path("${lib}.clust1.samples")

  script:
  """
  singularity run ${image} "freemuxlet \
    --plp ${plp}/${lib}.pileup \
    --nsample 18 \
    --out ${lib}"
  gunzip "${lib}.clust1.samples.gz"
  """
}

workflow {
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

  asis_ch = GenerateASISVCF(lib_ch, file(params.lines), driver_ch.dgrp_asis)

  // Mask highly expressed genes from informative variants
  mask_ch = MaskHighlyExpGenes(asis_ch, file(params.libpaths), file(params.gtf))

  // Permissively (> 300 fragments & > 200 genes) filter each library
  bclist_ch = PermissiveFilterBarcodes(lib_ch, file(params.libpaths))
  demux_vch = bclist_ch.join(mask_ch)

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

  // Calculate coverage per bam to remove SNPs without sufficient support
  coverage_ch = fbam_ch.combine(mask_ch, by: 0)
  cinfo_ch = FilterByCoverage(coverage_ch)

  mbam_ch = MergeBAM(fbam_ch)
  sinfo_ch = SortVCF(cinfo_ch.join(mbam_ch))

  // Run mpileup per 100 barcode to run in parallel
  sbclist_ch = bclist_ch
    | map { lib, bcpath ->
      def chunk = bcpath.splitText(by: 100, file: true)
      return [lib, chunk]
    }
    | transpose
  plp_ch = Mpileup(sinfo_ch.join(mbam_ch).combine(sbclist_ch, by: 0), file(params.popscle))

  intm_ch = plp_ch
    | groupTuple()
    | map { l, p ->
      def fp = p.flatten()
      return [l, fp]
    }
    | IntegratePL
  vplp_ch = intm_ch
    | map { l, p ->
      def pbase = p[0].getParent()
      return [l, pbase]
    }

  demux_out = Demuxlet(vplp_ch.join(cinfo_ch), params.popscle)
  freemux_out = Freemux(vplp_ch, params.popscle)

  demux_out | view
  freemux_out | view
}
