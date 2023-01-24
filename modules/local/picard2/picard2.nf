#!/usr/bin/env nextflow

process picard_create_seq_dict {
// Runs CreateSequenceDictionary
//
// input:
//   path(fa) - Reference FASTA
//   val(parstr) - Additional Parameters
//
// output:
//   tuple => emit: dict_file
//     path(fa) - Reference FASTA
//     path('*.dict') - Dict File
//
// require:
//   FA

  tag "${fa}"
  label 'picard2_container'
  label 'picard2_create_seq_dict'

  input:
  path fa
  val parstr

  output:
  tuple path(fa), path('*.dict'), emit: dict_file

  script:
  """
  fa_buf=`echo ${fa}`
  ofa=\${fa_buf%.fa*}.dict
  java -XX:ParallelGCThreads=${task.cpus} -jar /usr/picard/picard.jar CreateSequenceDictionary R=${fa} O=\${ofa} ${parstr}
  """
}


process picard_mark_duplicates {
// Runs MarkDuplicates
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: mkdup_bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*.mkdup.bam") - MarkDup Alignment File
//
// require:
//   BAMS
//   params.picard_mark_duplicates_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'picard2_container'
  label 'picard2_mark_duplicates'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/picard2_mark_duplicates", pattern: "*.marked_dup_metrics.txt"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.mkdup.bam"), emit: mkdup_bams
  tuple val(pat_name), val(run), val(dataset), path("*.marked_dup_metrics.txt"), emit: marked_dup_metrics

  script:
  """
  mkdir -p tmp
  bam_buf=`echo ${bam}`
  ofa=\${bam_buf%.bam}.mkdup.bam
  java -XX:ParallelGCThreads=${task.cpus} -Xmx${task.memory.toGiga()}g -jar /usr/picard/picard.jar MarkDuplicates I=${bam} O=\${ofa} M=${dataset}-${pat_name}-${run}.marked_dup_metrics.txt ${parstr} TMP_DIR=\${PWD}/tmp
  rm -rf tmp/
  """
}


process picard_collect_insert_size_metrics {
// Runs CollectInsertSizeMetrics
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - Alignment File
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: insert_metrics
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - Insert Size Metrics File
//
// require:
//   ALNS
//   params.picard_collect_insert_size_metrics_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'picard2_container'
  label 'picard2_collect_insert_size_metrics'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/picard2_collect_insert_size_metrics"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.txt'), emit: insert_size_metrics
  tuple val(pat_name), val(run), val(dataset), path('*.pdf'), emit: insert_size_metrics_pdfs

  script:
  """
  java -jar /usr/picard/picard.jar CollectInsertSizeMetrics \
    I=${bam} \
    O=${dataset}-${pat_name}-${run}-insert_size_metrics.txt \
    H=${dataset}-${pat_name}-${run}-insert_size_metrics.pdf \
    ${parstr}
  """
}



process picard_collect_rna_seq_metrics {
// Runs CollectRnaSeqMetrics
// input:
// output:
//
// require:
//   ALNS
//   params.picard2$picard_collect_rna_seq_metrics_parameters
//   params.picard2$picard_collect_rna_seq_metrics_ref_flat
//   params.picard2$picard_collect_rna_seq_metrics_rrna_ref
//   params.picard2$picard_collect_rna_seq_metrics_rna_strand_specificity

  tag "${dataset}/${pat_name}/${run}"
  label 'picard2_container'
  label 'picard2_collect_rna_seq_metrics'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/picard2_collect_rna_seq_metrics", pattern: "*rnaqc*"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  val parstr
  path ref_flat
  path rrna_ref
  val strand_specificity

  output:
  tuple val(pat_name), val(run), val(dataset), path(bam), path("*.txt"), emit: rna_seq_metrics_reports
  tuple val(pat_name), val(run), val(dataset), path(bam), path("*.pdf"), emit: rna_seq_metrics_report_pdfs

  script:
  """
  ABAM=`echo ${bam}`
  java -Xmx${task.memory.toGiga()}g \
    -XX:ParallelGCThreads=${task.cpus} \
    -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
    REF_FLAT=${ref_flat} \
    RIBOSOMAL_INTERVALS=${rrna_ref} \
    STRAND_SPECIFICITY=${strand_specificity} \
    CHART_OUTPUT=\${ABAM%.bam*}.rnaqc.pdf INPUT=${bam} \
    OUTPUT=\${ABAM%.bam*}.rnaqc.txt
  """
}


process picard_collect_wgs_metrics_nzc {
// Collect metrics about coverage and performance of whole genome sequencing
// (WGS) experiments. This tool collects metrics about the percentages of reads
// that pass base- and mapping- quality filters as well as coverage
// (read-depth) levels. Both minimum base- and mapping-quality values as well
// as the maximum read depths (coverage cap) are user defined. This extends
// CollectWgsMetrics by including metrics related only to siteswith non-zero
// (>0) coverage.
//
// input:
// output:
//
// require:
//   BAMS
//   params.picard2$picard_collect_wgs_metrics_nzc$dna_ref
//   params.picard2$picard_collect_wgs_metrics_nzc$collect_wgs_metrics_nzc_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'picard2_container'
  label 'picard2_collect_wgs_metrics_nzc'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/picard2_collect_wgs_metrics_nzc"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  path fa
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*metrics.txt"), emit: wgs_metrics
  tuple val(pat_name), val(run), val(dataset), path("*metrics.pdf"), emit: wgs_metrics_pdfs

  script:
  """
  java -XX:ParallelGCThreads=${task.cpus} -jar /usr/picard/picard.jar CollectWgsMetricsWithNonZeroCoverage I=${aln} O=${dataset}-${pat_name}-${run}.wgs_metrics.txt CHART=${dataset}-${pat_name}-${run}.wgs_metrics.pdf R=${fa}
  """
}


process picard_collect_vcf_metrics {
// Collects per-sample and aggregate (spanning all samples) metrics from the
// provided VCF file.
//
// input:
// output:
//
// require:
//   VCFS
//   params.picard2$picard_collect_vcf_metrics$dbsnp_ref
//   params.picard2$picard_collect_vcf_metrics$collect_vcf_metrics_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'picard2_container'
  label 'picard2_collect_vcf_metrics'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/picard2_collect_vcf_metrics"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf), path(csi)
  path dbsnp_ref
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*detail_metrics"), emit: detail_metrics
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*summary_metrics"), emit: summary_metrics

  script:
  """
  AVCF=`echo ${vcf}`
  java -XX:ParallelGCThreads=${task.cpus} -jar /usr/picard/picard.jar CollectVariantCallingMetrics I=${vcf} O=\${AVCF%.vcf*}.vcf_metrics.txt DBSNP=${dbsnp_ref}
  """
}
