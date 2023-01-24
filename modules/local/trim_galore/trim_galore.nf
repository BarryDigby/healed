#!/usr/bin/env nextflow

process trim_galore {
// Runs trim galore
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
  label 'trim_galore_container'
  label 'trim_galore'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/trim_galore"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*_1*.trimmed.f*q.gz"), path("${dataset}-${pat_name}-${run}*_2*.trimmed.f*q.gz"), emit: procd_fqs
  tuple val(pat_name), val(run), val(dataset), path("*_1_fastqc.zip"), path("*_2_fastqc.zip"), optional: true, emit: fastqc_zips

  script:
  """
  trim_galore --basename ${run} ${parstr} ${fq1} ${fq2} -j ${task.cpus}
  mv ${run}_R1_val_1.fq.gz ${dataset}-${pat_name}-${run}_1.trimmed.fq.gz
  mv ${run}_R2_val_2.fq.gz ${dataset}-${pat_name}-${run}_2.trimmed.fq.gz
  """
}

process trim_galore_hlap {
// Runs trim galore with hard trimming for HLAProfiler. Needs to be wrapped up
// into the primary trim_galore process definition.
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
  label 'trim_galore_container'
  label 'trim_galore'
//  publishDir "${params.samps_out_dir}/${dataset}/${run}/trim_galore"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}*_1*.trimmed*50bp_5prime.f*q.gz"), path("${dataset}-${pat_name}-${run}*_2*.trimmed*50bp_5prime.f*q.gz"), emit: procd_fqs

  script:
  """
  trim_galore --basename ${run} ${parstr} ${fq1} ${fq2} -j ${task.cpus}
  """
}
