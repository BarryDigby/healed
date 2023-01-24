#!/usr/bin/env nextflow

process fastqc {
// A quality control tool for high throughput sequence data.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - Fastq 1
//     path(fq2) - Fastq 2
//   val(parstr) - Additional Parameters
//
// output:
//   tuple => emit: fastqc_reports
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*") - Outputs
//   tuple => emit: fastqc_zips
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*zip") - Output Zip Files
//   tuple => emit: fastqc_htmls
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*htmls") - Output HTML Files
//
// required:
//   FQS
//   params.fastqc$fastqc_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'fastqc_container'
  label 'fastqc'
  publishDir "${params.qc_out_dir}/${dataset}/${pat_name}/${run}/fastqc"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*"), emit: fastqc_reports
  tuple val(pat_name), val(run), val(dataset), path("*zip"), emit: fastqc_zips
  tuple val(pat_name), val(run), val(dataset), path("*html"), emit: fastqc_htmls

  script:
  """
  fastqc ${parstr} ${fq1} ${fq2} -t ${task.cpus}
  """
}
