#!/usr/bin/env nextflow

process trimmomatic {
// Runs trimmomatic
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - FASTQ Prefix
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   val parstr - Parameter String
//
// output:
//   tuple => emit: procd_fqs
//     val(pat_name) - Patient Name
//     val(run) - FASTQ Prefix
//     val(dataset) - Dataset
//     path("${dataset}-${pat_name}-${run}*_1*.trimmed.f*q.gz") - Trimmed FASTQ 1
//     path("${dataset}-${pat_name}-${run}*_2*.trimmed.f*q.gz") - Trimmed FASTQ 2
//   path("meta") - Metadata File

// require:
//   FQS
//   params.trim_galore$trim_galore_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'trimmomatic_container'
  label 'trimmomatic'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/trimmomatic"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*_1*.trimmed.f*q.gz"), path("${dataset}-${pat_name}-${run}*_2*.trimmed.f*q.gz"), emit: procd_fqs
  tuple val(pat_name), val(run), val(dataset), path("*_1_fastqc.zip"), path("*_2_fastqc.zip"), optional: true, emit: fastqc_zips

  script:
  def SE_OR_PE =  fq2 == '' ? 'SE' : 'PE'
  def RD2 = fq2 != '' ? "${dataset}-${pat_name}-${run}_2.trimmed.fq.gz ${dataset}-${pat_name}-${run}_2.orphan.fq.gz" : ""
  """
  trimmomatic $SE_OR_PE -threads ${task.cpus} ${fq1} ${fq2} \
  ${dataset}-${pat_name}-${run}_1.trimmed.fq.gz \
  ${dataset}-${pat_name}-${run}_1.orphan.fq.gz \
  $RD2 \
  ${parstr}
  """
}
