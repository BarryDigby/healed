#!/usr/bin/env nextflow


process pyclonevi_fit{
// Runs pyclone-vi fit
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq1) - FASTQ 2
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//   parstr - Additional Parameters
//
// output:
//   tuple => emit: sams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.sam') - Alignment Output File

// require:
//   PCVIS
//   params.pyclone-vi$pyclone-vi_fit_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'pyclonevi_container'
  label 'pyclonevi_fit'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/pyclone-vi_fit"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset),  path(pcvi)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*.tmp'), optional: true, emit: pcvi_tmps

  script:
  """
  source activate pyclone-vi
  pyclone-vi fit -i ${pcvi} -o ${pcvi}.tmp
  """
}

process pyclonevi_write_results_file {
// Runs pyclone-vi write_results_file
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq1) - FASTQ 2
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//   parstr - Additional Parameters
//
// output:
//   tuple => emit: sams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path('*.sam') - Alignment Output File

// require:
//   PCVIS
//   params.pyclone-vi$pyclone-vi_write_results_file_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'pyclonevi_container'
  label 'pyclonevi_write_results_file'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/pyclone-vi_write_results_file"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset),  path(pcvi_tmp)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path('*.results'), emit: pcvi_results

  script:
  """
  source activate pyclone-vi
  APCVI=`echo ${pcvi_tmp}`
  pyclone-vi write-results-file -i ${pcvi_tmp} -o \${APCVI%.tmp}.results
  """
}
