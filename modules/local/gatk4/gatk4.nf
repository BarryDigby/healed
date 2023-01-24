#!/usr/bin/env nextflow

include { extract_chroms_from_bed } from '../ref_utils/ref_utils.nf'
include { bcftools_index_somatic } from '../bcftools/bcftools.nf'
include { bcftools_concat } from '../bcftools/bcftools.nf'
include { htslib_bgzip_somatic } from '../htslib/htslib.nf'

process gatk_haplotypecaller {
// Runs gatk HaplotypeCaller on a sample
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     path(bam) - BAM File
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//     path(dict_file) - Sequence Dictionary File
//   parstr - Additional Parameters
//   suffix - Output File Suffix
//
// output:
//   tuple => emit: germline_vcfs
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*germline.vcf") - Variant Call File

// require:
//   BAMS
//   REFS_W_INDICES
//   params.gatk$haplotypecaller$haplotypecaller_parameters
//   params.gatk$haplotypecaller$haplotypecaller_suffix

  tag "${dataset}/${pat_name}/${run}"
  label 'gatk4_container'
  label 'gatk4_haplotypecaller'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/gatk4_haplotypecaller"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  tuple path(fa), path(idx_files), path(dict_file)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(run), val(dataset), path("*germline.vcf"), emit: germline_vcfs

  script:
  """
  gatk HaplotypeCaller \
    -R ${fa} \
    -I ${bam} \
    ${parstr} \
    -O ${dataset}-${pat_name}-${run}.${suffix}.germline.vcf
  """
}


process gatk_mutect2_matched_by_chr {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'gatk4_container'
  label 'gatk4_mutect'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/gatk4_mutect2"

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai), val(chr)
  tuple path(fa), path(faidx_file), path(dict_file)
  path bed
  val parstr
  val suffix
  tuple path(pon_vcf), path(pon_vcf_tbi)
  tuple path(af_vcf), path(af_vcf_tbi)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.vcf"), emit: vcfs
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.vcf"), path("*stats*"), emit: vcfs_w_stats
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.f1r2.tar.gz"), emit: f1r2_tar_gzs

  script:
  def pon_proxy  = pon_vcf.name != 'dummy_file' ? "-pon ${pon_vcf}" : ''
  def af_proxy  = af_vcf.name != 'dummy_file' ? "--germline-resource ${af_vcf}" : ''


  """
  CHR=`echo ${chr}`
  grep "\${CHR%\n}	" ${bed} > \${CHR}.bed

  gatk Mutect2 \
    -L \${CHR}.bed \
    -R ${fa} \
    -I ${tumor_bam} \
    -I ${norm_bam} \
    -tumor ${dataset}-${pat_name}-${tumor_run} \
    -normal ${dataset}-${pat_name}-${norm_run} \
    ${af_proxy} \
    ${pon_proxy} \
    --f1r2-tar-gz ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${CHR}.f1r2.tar.gz \
    ${parstr} \
    -O ${dataset}-${pat_name}-${norm_run}_${tumor_run}.\${CHR}${suffix}.vcf
  """
}


workflow gatk_mutect2_matched {
  take:
    bams
    ref
    bed
    parstr
    suffix
    pon_vcf
    af_vcf
    species
  main:
    if( species =~ /[Hh]uman|hs|HS|[Hh]omo/ ) {
        chrom_count = 25
    } else if( species =~ /[Mm]ouse|mm|MM|[Mm]us/ ) {
        chrom_count = 21
    }
    extract_chroms_from_bed(
      bed)
    bams
      .combine(extract_chroms_from_bed.out.chroms_list.splitText())
      .set{ bams_w_chr }
    gatk_mutect2_matched_by_chr(
      bams_w_chr,
      ref,
      bed,
      parstr,
      suffix,
      pon_vcf,
      af_vcf)
    htslib_bgzip_somatic(
      gatk_mutect2_matched_by_chr.out.vcfs)
    bcftools_index_somatic(
      htslib_bgzip_somatic.out.bgzip_files,
      '')
    bcftools_index_somatic.out.vcfs_w_csis
      .groupTuple(by: [0, 1, 2, 3], size: chrom_count)
      .set{ vcfs_by_patient }
    bcftools_concat(
      vcfs_by_patient,
      '',
      'mutect2')
    gatk_mutect2_matched_by_chr.out.vcfs_w_stats
      .map{ [it[0], it[1], it[2], it[3], it[5]] }
      .groupTuple(by: [0, 1, 2, 3], size: chrom_count)
      .set{ vcf_stats_by_patient }
    gatk_merge_mutect_stats(
      vcf_stats_by_patient,
      '')
    bcftools_concat.out.concatd_vcfs
      .join(gatk_merge_mutect_stats.out.merged_vcf_stats, by: [0, 1, 2, 3])
      .set { vcfs_w_stats }
  emit:
    vcfs = bcftools_concat.out.concatd_vcfs
    vcfs_w_stats = vcfs_w_stats
    f1r2_tar_gzs = gatk_mutect2_matched_by_chr.out.f1r2_tar_gzs
    
}


process gatk_base_recalibrator {
// Runs gatk BaseRecalibrator
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File
//     path(bai) - BAI File
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//     path(dict_file) - Sequence Dictionary File
//   tuple
//     path(known_sites) - Known Sites Reference VCF
//     path(known_sites_index) - Known Sites Reference VCF Index
//   path bed - BED File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: grps
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File (same as input)
//     path("*grp") - Recalibration Matrix

// require:
//   ALNS
//   REF_W_INDICES
//   KNOWN_SITES_W_INDEX
//   BED
//   params.gatk$base_recalibrator_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'gatk4_container'
  label 'gatk4_base_recalibrator'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  tuple path(fa), path(idx_files), path(dict_file)
  tuple path(known_sites), path(known_sites_index)
  path bed
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path("*grp"), emit: grps

  script:
  """
  gatk BaseRecalibrator \
    -I ${bam} \
    -R ${fa} \
    -L ${bed} \
    --known-sites ${known_sites} \
    -O ${dataset}-${pat_name}-${run}.recal_data.grp ${parstr}
  """
}


process gatk_apply_bqsr {
// Runs gatk ApplyBQSR
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bed) - BED File
//     path(grp) - Recalibration Matrix
//   tuple
//     path(fa) - Reference FASTA
//     //idx1 and idx2 can be removed and the incoming channel can be truncated
//     //before being passed to gatk_apply_bqsr.
//     path(idx1) - ???
//     path(idx2) - ???
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*bam") - Output BAM File

// require:
//   ALNS_W_GRPS
//   REF_W_INDICES
//   params.gatk$apply_bqsr_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'gatk4_container'
  label 'gatk4_apply_bqsr'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai), path(grp)
  tuple path(fa), path(idx1), path(idx2)
  path bed
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*bam"), emit: bams

  script:
  """
  ALN_NAME=`echo ${bam}`
  gatk ApplyBQSR \
    -L ${bed} \
    -I ${bam} \
    -R ${fa} \
    --bqsr-recal-file ${grp} \
    -O \${ALN_NAME%.bam}.bqsr.bam ${parstr}
  """
}


process gatk_index_feature_file {
// Runs gatk IndexFeatureFile
//
// input:
//   path feature_file - Feature File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: ff_w_index
//     path(feature_file) - Feature File
//     path("${feature_file}.tbi"}) - Feature File Index

// require:
//   params.gatk_index_feature_file$feature_file
//   params.gatk_index_feature_file$gatk_index_feature_file_parameters

  tag "${feature_file}"
  label 'gatk4_container'
  label 'gatk4_index_feature_file'

  input:
  path feature_file
  val parstr

  output:
  tuple path(feature_file), path("${feature_file}.tbi"), emit: ff_w_index

  script:
  def ff_proxy  = feature_file.name != 'dummy_file' ? "gatk IndexFeatureFile ${parstr} -I ${feature_file}" : "touch ${feature_file}.tbi"
  """
  ${ff_proxy}
  """
}


process gatk_merge_mutect_stats {
// Runs gatk MergeMutectStats on VCF stats files.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(dataset) - Dataset
//     val(norm_run) - Normal Run Name
//     path(norm_bam) - Normal BAM
//     path(norm_bai) - Normal BAI
//     val(tumor_run) - Tumor Run Name
//     path(tumor_bam) - Tumor BAM
//     path(tumor_bai) - Tumor BAI
//   tuple
//     path(fa) - Reference FASTA
//     path(faidx_file) - FAIDX File
//     path(dict_file) - Sequence Dictionary File
//   val parstr - Additional Paramters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*.vcf")

// require:
//   VCF_STATS
//   params.gatk$gatk_filter_mutect_calls_parameters
//   params.gatk$gatk_filter_mutect_calls_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'gatk4_container'
  label 'gatk4_merge_mutect_stats'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/gatk4_merge_mutect_stats"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf_stats)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.merged.vcf.stats"), emit: merged_vcf_stats

  script:
  """
  STATS=\$(ls -v *vcf.stats)
  STATS_STRING=''
  for stats in \${STATS}; do
    STATS_STRING="\${STATS_STRING} -stats \${stats}"
  done
  gatk MergeMutectStats \${STATS_STRING} -O ${dataset}-${pat_name}-${norm_run}_${tumor_run}.merged.vcf.stats
  """
}


process gatk_filter_mutect_calls {
// Runs gatk FilterMutectCalls on VCFs.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(dataset) - Dataset
//     val(norm_run) - Normal Run Name
//     path(norm_bam) - Normal BAM
//     path(norm_bai) - Normal BAI
//     val(tumor_run) - Tumor Run Name
//     path(tumor_bam) - Tumor BAM
//     path(tumor_bai) - Tumor BAI
//   tuple
//     path(fa) - Reference FASTA
//     path(faidx_file) - FAIDX File
//     path(dict_file) - Sequence Dictionary File
//   val parstr - Additional Paramters
//   val suffix - Output File Suffix
//
// output:
//   tuple => emit: vcfs
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*.vcf")

// require:
//   ALNS
//   REF_W_INDICES
//   params.gatk$gatk_filter_mutect_calls_parameters
//   params.gatk$gatk_filter_mutect_calls_suffix

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'gatk4_container'
  label 'gatk4_filter_mutect_calls'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/gatk4_filter_mutect_calls"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(vcf_stats), path(contamination_table), path(artifact_priors)
  tuple path(fa), path(faidx_file), path(dict_file)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.vcf"), emit: filtd_vcfs

  script:
  def contam_proxy = contamination_table.name != 'dummy_file' ? "--contamination-table ${contamination_table}" : ""
  """
  ORIG="${vcf}"
  gatk FilterMutectCalls \
    -R ${fa} \
    -V ${vcf} \
    --stats ${vcf_stats} \
    ${contam_proxy} \
    --orientation-bias-artifact-priors ${artifact_priors} \
    -O \${ORIG%.vcf*}${suffix}.vcf
  """
}


process gatk_get_pileup_summaries {
// Runs gatk GetPileupSummaries on BAMs
//
  
  tag "${dataset}/${pat_name}/${run}"
  label 'gatk4_container'
  label 'gatk4_get_pileup_summaries'
  
  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  tuple path(sites_file), path(sites_file_idx)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}.pileups.table"), emit: pileups_tables

  script:
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}G" GetPileupSummaries \
    -I ${bam} \
    -V ${sites_file} \
    -L ${sites_file} \
    ${parstr} \
    -O ${dataset}-${pat_name}-${run}.pileups.table
  """
}


process gatk_calculate_contamination {
// Runs gatk CalculateContamination
  
  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'gatk4_container'
  label 'gatk4_calculate_contamiantion'
  
  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_pileups_table), val(tumor_run), path(tumor_pileups_table)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("${dataset}-${pat_name}-${norm_run}_${tumor_run}.contamination_table"), emit: contamination_tables

  script:
  """
  gatk CalculateContamination \
    -I ${tumor_pileups_table} \
    -matched ${norm_pileups_table} \
    -tumor-segmentation ${dataset}-${pat_name}-${norm_run}_${tumor_run}.segments.table \
    -O ${dataset}-${pat_name}-${norm_run}_${tumor_run}.contamination_table
  """
}


process gatk_learn_read_orientation_model {
// Runs LearnReadOrientationModel
//
  
  tag "${dataset}/${pat_name}/${run}"
  label 'gatk4_container'
  label 'gatk4_learn_read_orientation_model'
  
  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(all_f1r2s)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("${dataset}-${pat_name}-${norm_run}_${tumor_run}.artifact_prior.tar.gz"), emit: artifact_priors

  script:
  """
  f1r2s_STR=""
  for i in `ls *f1r2.tar.gz`; do
    f1r2s_STR="\${f1r2s_STR}-I \${i} "
  done

  gatk LearnReadOrientationModel \${f1r2s_STR} -O ${dataset}-${pat_name}-${norm_run}_${tumor_run}.artifact_prior.tar.gz
  """
}
