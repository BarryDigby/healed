#!/usr/bin/env nextflow

//include { qc_rna_seq } from '../binfotron/binfotron.nf'
//include { qc_fastq } from '../binfotron/binfotron.nf'
include { fastqc } from '../fastqc/fastqc.nf'
include { kraken2 } from '../kraken2/kraken2.nf'
include { multiqc as multiqc_subd } from '../multiqc/multiqc.nf'
include { multiqc as multiqc_subd_by_pat } from '../multiqc/multiqc.nf'
include { seqtk_sample as seqtk_sample_dna } from '../seqtk/seqtk.nf'
include { seqtk_sample as seqtk_sample_rna } from '../seqtk/seqtk.nf'
include { trim_galore as trim_galore_dna } from '../trim_galore/trim_galore.nf'
include { trim_galore as trim_galore_rna } from '../trim_galore/trim_galore.nf'
include { manifest_to_raw_fqs as manifest_to_raw_fqs_dna } from '../preproc/preproc.nf'
include { manifest_to_raw_fqs as manifest_to_raw_fqs_rna } from '../preproc/preproc.nf'
include { picard_collect_rna_seq_metrics } from '../picard2/picard2.nf'
include { picard_collect_wgs_metrics_nzc } from '../picard2/picard2.nf'
include { picard_collect_vcf_metrics } from '../picard2/picard2.nf'
include { procd_fqs_to_raw_alns } from '../somatic/somatic.nf'
include { gtf_to_refflat } from '../ref_utils/ref_utils.nf'
include { gtf_to_rrna_intervals } from '../ref_utils/ref_utils.nf'
include { procd_fqs_to_star_alns } from '../rna_quant/rna_quant.nf'
include { samtools_sort } from '../samtools/samtools.nf'
include { get_fastqs } from '../utilities/utilities.nf'
include { samtools_stats } from '../samtools/samtools.nf'

include { lenstools_consolidate_multiqc_stats } from '../lenstools/lenstools.nf'
//include { binfo_qc } from '../binfotron/binfotron.nf'


workflow bams_to_rna_seq_metrics {
// require:
//   params.qc$bams_to_rna_seq_metrics$qc_genome
//   BAMS
//   params.qc$bams_to_rna_seq_metrics$gtf
//   params.qc$bams_to_rna_seq_metrics$picard_collect_rna_seq_metrics_parameters
//   params.qc$bams_to_rna_seq_metrics$picard_strand_specificity
  take:
    fa
    bams
    gtf
    picard_collect_rna_seq_metrics_parameters
    picard_strand_specificity
    samtools_faidx_parameters
  main:
    gtf_to_refflat(
      gtf)
    gtf_to_rrna_intervals(
      fa,
      gtf,
      samtools_faidx_parameters)
    picard_collect_rna_seq_metrics(
      bams,
      picard_collect_rna_seq_metrics_parameters,
      gtf_to_refflat.out.refflat,
      gtf_to_rrna_intervals.out.rrna_intervals,
      picard_strand_specificity)
  emit:
    rna_seq_metrics_reports = picard_collect_rna_seq_metrics.out.rna_seq_metrics_reports
    //rna_seq_metrics_pdfs = picard_collect_rna_seq_metrics.out.rna_seq_metrics_pdfs
}


workflow fastqs_to_fastqc {
// require:
//   FASTQS
  take:
    fastqs
  main:
    fastqc(
      fastqs,
      params.qc$fastqs_to_fastqc$fastqc_parameters)
  emit:
    fastqc_reports = fastqc.out.fastqc_reports
}


workflow manifest_to_kraken2 {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest.map{ [it[0], it[1], it[2], it[3]] },
      params.fq_dir)
    raw_fqs_to_kraken2(
      get_fastqs.out.fastqs)
  emit:
    kraken_outs = raw_fqs_to_kraken2.out.kraken_outs
}


workflow procd_fqs_to_kraken2 {
// require:
//   fqs
  take:
    fqs
  main:
    kraken2(
      fqs,
      params.kraken2$procd_fqs_to_kraken2$kraken2_database,
      params.kraken2$procd_fqs_to_kraken2$kraken2_parameters)
  emit:
    kraken_outs = kraken2.out.kraken_outs
}


workflow filter_samples {
// require:
//   params.qc$filter_samp_ids$samps
//   params.qc$filter_samp_ids$qc_pass_files
  take:
    samps
    qc_pass_files
  main:
    filter_samples_sub(samps, qc_pass_files)
    filter_samples_sub.out.passes.map{ it -> [it[0], it[1], it[2], it[3]]}.set{ filtered_samps }
  emit:
    filtered_samps = filtered_samps
}


process filter_samples_sub {
// Strategy here is to filter samples by ensuring they are contained within all
// QC files, creating a flag file, and then not emitting that sample if the
// flag file doesn't exist.

  input:
  tuple val(pat_name), val(run), val(dataset), path(fle)
  path qc_pass_files

  output:
  tuple val(pat_name), val(run), val(dataset), path(fle), path('flag.file'), emit: passes optional true

  script:
  """
  QC_FILES_ACTUAL_COUNT=\$(ls *samples | wc -l)
  QC_FILES_DETECTED_COUNT=\$(grep ${run} *samples | wc -l)

  if [ \${QC_FILES_ACTUAL_COUNT} -eq \${QC_FILES_DETECTED_COUNT} ]; then
    touch 'flag.file'
  fi
  """
}


workflow bams_to_wgs_nzc_metrics {
// require:
//   ALNS
//   params.qc$bams_to_wgs_nzc_metrics$dna_ref
  take:
    alns
    fa
  main:
    picard_collect_wgs_metrics_nzc(
      alns,
      fa,
      params.qc$bams_to_wgs_nzc_metrics$collect_wgs_nzc_parameters)
//  emit:
}


workflow vcfs_to_vcf_metrics {
// require:
//   VCFS
//   params.qc$vcfs_to_vcf_metrics$dbsnp_ref
  take:
    vcfs
    dbsnp
    picard_collect_vcf_metrics_parameters
  main:
    picard_collect_vcf_metrics(
      vcfs,
      dbsnp,
      picard_collect_vcf_metrics_parameters)
  emit:
    detail_metrics = picard_collect_vcf_metrics.out.detail_metrics
    summary_metrics = picard_collect_vcf_metrics.out.summary_metrics
}


workflow initial_qc {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest.map{ [it[0], it[1], it[2], it[3]] },
      params.fq_dir)

    // DNA QC 
    manifest_to_raw_fqs_dna(
      manifest.filter{ it[4] =~ /WES/ })
    seqtk_sample_dna(
      manifest_to_raw_fqs_dna.out.fastqs,
      params.qc$filter_samps_by_qc$seqtk_sample_count,
      params.qc$filter_samps_by_qc$seqtk_sample_seed,
      params.qc$filter_samps_by_qc$seqtk_sample_suffix,
      params.qc$filter_samps_by_qc$seqtk_sample_parameters)
    trim_galore_dna(
      seqtk_sample_dna.out.subd_fqs,
      params.qc$filter_samps_by_qc$trim_galore_parameters)
    procd_fqs_to_raw_alns(
      params.qc$filter_samps_by_qc$dna_ref,
      trim_galore_dna.out.procd_fqs,
      manifest)
    samtools_stats(
      procd_fqs_to_raw_alns.out.alns,
      '')
    

    // RNA QC
    manifest_to_raw_fqs_rna(
      manifest.filter{ it[4] =~ /RNA/ })
    seqtk_sample_rna(
      manifest_to_raw_fqs_rna.out.fastqs,
      params.qc$filter_samps_by_qc$seqtk_sample_count,
      params.qc$filter_samps_by_qc$seqtk_sample_seed,
      params.qc$filter_samps_by_qc$seqtk_sample_suffix,
      params.qc$filter_samps_by_qc$seqtk_sample_parameters)
    trim_galore_rna(
      seqtk_sample_rna.out.subd_fqs,
      params.qc$filter_samps_by_qc$trim_galore_parameters)
    procd_fqs_to_star_alns(
      params.qc$filter_samps_by_qc$star_ref,
      trim_galore_rna.out.procd_fqs)
//    bams_to_rna_seq_metrics(
//      params.qc$filter_samps_by_qc$star_ref,
//      procd_fqs_to_star_alns.out.alns,
//      params.qc$filter_samps_by_qc$gtf)



//    trim_galore_dna.out.fastqc_zips.concat(trim_galore_rna.out.fastqc_zips).concat(samtools_stats.out.bam_stats).concat(bams_to_rna_seq_metrics.out.rna_seq_metrics_reports).set{ all_metrics_files }
    trim_galore_dna.out.fastqc_zips.concat(trim_galore_rna.out.fastqc_zips).concat(samtools_stats.out.bam_stats).set{ all_metrics_files }

//    fastqs_to_fastqc.out.fastqc_reports.concat(samtools_stats.out.bam_stats).concat(bams_to_rna_seq_metrics.out.rna_seq_metrics_reports).set{ all_metrics_files }
//    all_metrics_files.groupTuple(by: 0).map{ it -> it[3] }.set{ metrics_files_by_pat }

    multiqc_subd(
     all_metrics_files.map{ it -> it[3] }.collect(),
      '')

//    multiqc_subd_by_pat(
//     metrics_files_by_pat,
//      '')
//  emit

  lenstools_consolidate_multiqc_stats(
    multiqc_subd.out.multiqc_data)
}
