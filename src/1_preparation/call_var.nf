params.fasta = '/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.fa'
params.dict = '/scratch/cgsb/desplan/File_exchange/Yen_ref/ensembl_88_Nikos/Drosophila_melanogaster.BDGP6.dna_sm.toplevel.dict'
params.bamdir = './int/wgs_bam/*.bam'

process AddRG {
  cpus '1'
  memory '4GB'
  time '2h'
  module 'picard/2.27.5'

  input:
  path bam

  output:
  path "${bam.baseName}_RG.bam"

  script:
  """
  java -Xmx3G -jar \${PICARD_JAR} AddOrReplaceReadGroups \
    I=${bam} \
    O=${bam.baseName}_RG.bam \
    SORT_ORDER=coordinate \
    RGLB=lib1 \
    RGPL=illumina \
    RGPU=unit1 \
    RGSM=${bam.baseName} \
    CREATE_INDEX=False
  """
}

process MarkDuplicate {
  cpus '1'
  memory '4GB'
  time '2h'
  module 'picard/2.27.5'

  input:
  path bam

  output:
  path "${bam.baseName}_markdup.bam"

  script:
  """
  java -Xmx3G -jar \${PICARD_JAR} MarkDuplicates \
    I=${bam} \
    O=${bam.baseName}_markdup.bam \
    M=marked_dup_metrics.txt
  """
}

process IndexBam {
  cpus '1'
  memory '4GB'
  time '1h'
  module 'picard/2.27.5'

  input:
  path bam

  output:
  tuple path("${bam}"), path("${bam.baseName}.bai")

  script:
  """
  java -Xmx4G -jar \${PICARD_JAR} BuildBamIndex \
    I=${bam}
  """
}

process CallHaplotype {
  cpus '4'
  memory '4GB'
  time '8h'
  module 'gatk/4.3.0.0'

  input:
  tuple path(bam), path(bai), val(chr)
  path ref
  path ref_dict
  path ref_index

  output:
  tuple path("${bam.baseName}_${chr}.gvcf.gz"), path("${bam.baseName}_${chr}.gvcf.gz.tbi")

  script:
  """
  gatk HaplotypeCaller \
    --java-options '-Xmx3G' \
    -R ${ref} \
    -I ${bam} \
    -L ${chr} \
    --native-pair-hmm-threads ${task.cpus} \
    -O "${bam.baseName}_${chr}.gvcf.gz" \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    -G AS_StandardAnnotation \
    -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90 \
    -ERC GVCF
  """
}

process CombineGVCFs {
  cpus '1'
  memory '4GB'
  time '1h'
  module 'gatk/4.3.0.0'
  publishDir 'int/wgs_variant/', mode: 'copy', overwrite: false

  input:
  path vcfs
  path tbis
  path ref
  path ref_dict
  path ref_index

  output:
  tuple path('merged.gvcf.gz'), path('merged.gvcf.gz.tbi')

  script:
  // Create a string of -V arguments for each GVCF file
  def variants = vcfs.collect { "-V ${it}" }.join(' ')
  """
  gatk CombineGVCFs \
    --java-options '-Xmx3G' \
    ${variants} \
    -R ${ref} \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    -G AS_StandardAnnotation \
    -O merged.gvcf.gz
  """
}

process JointGenotyping {
  cpus '1'
  memory '3GB'
  time '30m'
  module 'gatk/4.3.0.0'

  input:
  tuple path(gvcf), path(tbi)
  path ref
  path ref_dict
  path ref_index

  output:
  tuple path('joint_geno.vcf.gz'), path('joint_geno.vcf.gz.tbi')

  script:
  """
  gatk GenotypeGVCFs \
    --java-options "-Xmx2G" \
    -R ${ref} \
    -O 'joint_geno.vcf.gz' \
    -G StandardAnnotation \
    -G StandardHCAnnotation \
    -G AS_StandardAnnotation \
    -V ${gvcf}
  """
}

process INDELCalibrate {
  cpus '1'
  memory '2GB'
  time '1h'
  module 'gatk/4.3.0.0'

  input:
  tuple path(vcf), path(tbi)
  path ref
  path ref_dict
  path ref_index
  path refvar

  output:
  tuple path('VQSR_INDEL.recal'), path('VQSR_INDEL.tranches')

  script:
  def res = refvar[0]
  """
  gatk VariantRecalibrator \
    --java-options "-Xmx7G" \
    -R ${ref} \
    -V ${vcf} \
    -O VQSR_INDEL.recal \
    --tranches-file VQSR_INDEL.tranches \
    -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    --resource:hardfilter,known=true,training=true,truth=true,prior=15.0 ${res} \
    -mode INDEL
  """
}

process SNPCalibrate {
  cpus '1'
  memory '2GB'
  time '1h'
  module 'gatk/4.3.0.0'

  input:
  tuple path(vcf), path(tbi)
  path ref
  path ref_dict
  path ref_index
  path refvar

  output:
  tuple path('VQSR_SNP.recal'), path('VQSR_SNP.tranches')

  script:
  def res = refvar[0]
  """
  gatk VariantRecalibrator \
    --java-options "-Xmx7G" \
    -R ${ref} \
    -V ${vcf} \
    -O VQSR_SNP.recal \
    --tranches-file VQSR_SNP.tranches \
    -an DP -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    --resource:hardfilter,known=true,training=true,truth=true,prior=15.0 ${res} \
    -mode SNP
  """
}

process BcftoolsCall {
  cpus '16'
  memory '2GB'
  time '8h'
  module 'bcftools/intel/1.19'
  publishDir 'int/wgs_variant', mode: 'copy', overwrite: false

  input:
  path bam
  path bai
  path ref
  path ref_index

  output:
  path 'raw_variants.vcf.gz'

  script:
  def bams = bam.collect().join(' ')
  def cores = task.cpus / 2
  """
  bcftools mpileup \
    -f ${ref} \
    ${bams} \
    -d 8000 \
    --threads ${cores} |\
  bcftools call -mv -Oz \
    --threads ${cores} \
    -o 'raw_variants.vcf.gz'
  """
}

process GATKHardFilter {
  cpus '1'
  memory '2GB'
  time '1h'
  module 'gatk/4.3.0.0'

  input:
  tuple path(vcf), path(tbi)

  output:
  tuple path('gatk_snp.vcf.gz'), path('gatk_indel.vcf.gz')

  script:
  """
  gatk SelectVariants \
    --java-options "-Xmx3G" \
    -V ${vcf} \
    -select-type SNP \
    -O snps.vcf.gz

  gatk VariantFiltration \
    --java-options "-Xmx3G" \
    -V snps.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O gatk_snp.vcf.gz

  gatk SelectVariants \
    --java-options "-Xmx3G" \
    -V ${vcf} \
    -select-type INDEL \
    -O indels.vcf.gz

  gatk VariantFiltration \
    --java-options "-Xmx3G" \
    -V indels.vcf.gz \
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "FS > 200.0" --filter-name "FS200" \
    -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    -O gatk_indel.vcf.gz
  """
}
process MergeVCF {
  cpus '1'
  memory '3GB'
  time '1h'
  module 'picard/2.27.5'

  input:
  path vcfs

  output:
  path 'gatk_merged_filterd.vcf.gz'

  script:
  def vcf = vcfs.collect { v -> "I=${v}" }.join(' ')
  """
  java -Xmx3G -jar \${PICARD_JAR} MergeVcfs \
    ${vcf} \
    O=gatk_merged_filterd.vcf.gz
  """
}

process BcftoolsFilter {
  cpus '4'
  memory '2GB'
  time '1h'
  module 'bcftools/intel/1.19'

  input:
  path vcf

  output:
  path 'bcf_filtered.vcf.gz'

  script:
  """
  bcftools view \
    --threads ${task.cpus} \
    -i 'QUAL>=20' \
    -Oz \
    -o bcf_filtered.vcf.gz \
    ${vcf}
  """
}

process BcftoolsIntersect {
  cpus '1'
  memory '2GB'
  time '1h'
  module 'bcftools/intel/1.19'

  input:
  path vcf1
  path vcf2

  output:
  tuple path('output/0002.vcf.gz'), path('output/0002.vcf.gz.tbi')

  script:
  """
  bcftools index ${vcf1}
  bcftools index ${vcf2}
  bcftools isec -p output \
    -Oz ${vcf1} ${vcf2}
  """
}

process ApplyINDELVQSR {
  publishDir 'int/wgs_variant/', mode: 'copy', overwrite: false
  cpus '1'
  memory '2GB'
  time '30m'
  module 'gatk/4.3.0.0'

  input:
  tuple path(vcf), path(tbi)
  tuple path(recal), path(tranches)
  path ref
  path ref_dict
  path ref_index

  output:
  tuple path('gatk_indel.recalibrated.vcf.gz'), path('gatk_indel.recalibrated.vcf.gz.tbi')

  script:
  """
  gatk SelectVariants \
    -V ${vcf} \
    -select-type INDEL \
    -O indel.vcf.gz

  gatk IndexFeatureFile \
    -I indel.vcf.gz

  gatk IndexFeatureFile \
    -I ${recal}

  gatk ApplyVQSR \
    -R ${ref} \
    -V indel.vcf.gz \
    --recal-file ${recal} \
    --tranches-file ${tranches} \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode INDEL \
    -O gatk_indel.recalibrated.vcf.gz
  """
}

process ApplySNPVQSR {
  publishDir 'int/wgs_variant/', mode: 'copy', overwrite: false
  cpus '1'
  memory '2GB'
  time '30m'
  module 'gatk/4.3.0.0'

  input:
  tuple path(vcf), path(tbi)
  tuple path(recal), path(tranches)
  path ref
  path ref_dict
  path ref_index

  output:
  tuple path('gatk_snp.recalibrated.vcf.gz'), path('gatk_snp.recalibrated.vcf.gz.tbi')

  script:
  """
  gatk SelectVariants \
    -V ${vcf} \
    -select-type SNP \
    -O snps.vcf.gz

  gatk IndexFeatureFile \
    -I snps.vcf.gz

  gatk IndexFeatureFile \
    -I ${recal}

  gatk ApplyVQSR \
    -V snps.vcf.gz \
    -R ${ref} \
    --recal-file ${recal} \
    --tranches-file ${tranches} \
    --truth-sensitivity-filter-level 99.0 \
    --create-output-variant-index true \
    -mode SNP \
    -O gatk_snp.recalibrated.vcf.gz
  """
}
workflow {
  bam_ch = channel.fromPath(params.bamdir)
  chr_ch = channel.from('X', 'Y', '2R', '2L', '3R', '3L', '4')

  bam_ch
    | AddRG
    | MarkDuplicate
    | IndexBam
    | set { in_ch }

  in_chr_ch = in_ch.combine(chr_ch)

  gvcfi_ch = CallHaplotype(in_chr_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'))
  gvcf_ch = gvcfi_ch.multiMap { v ->
    gvcf: v[0]
    tbi: v[1]
  }

  // Calling with GATK
  gvcf_list = gvcf_ch.gvcf.toList()
  tbi_list = gvcf_ch.tbi.toList()
  mvcf_ch = CombineGVCFs(gvcf_list, tbi_list, file(params.fasta), file(params.dict), file(params.fasta + '.fai'))
  vcf_ch = JointGenotyping(mvcf_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'))

  // Calling with bcftools
  input_ch = in_ch.multiMap { v ->
    bam: v[0]
    bai: v[1]
  }

  rvcf_ch = BcftoolsCall(
    input_ch.bam.collect(),
    input_ch.bai.collect(),
    file(params.fasta),
    file(params.fasta + '.fai'),
  )

  vcf_ch
    | GATKHardFilter
    | MergeVCF
    | set { gatkhf_ch }

  rvcf_ch
    | BcftoolsFilter
    | set { bcft_ch }

  ref_ch = BcftoolsIntersect(gatkhf_ch, bcft_ch)

  snp_recal_ch = SNPCalibrate(vcf_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'), ref_ch)
  indel_recal_ch = INDELCalibrate(vcf_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'), ref_ch)
  ApplySNPVQSR(vcf_ch, snp_recal_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'))
  ApplyINDELVQSR(vcf_ch, indel_recal_ch, file(params.fasta), file(params.dict), file(params.fasta + '.fai'))
}
