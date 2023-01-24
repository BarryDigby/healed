#!/usr/bin/env nextflow


process htslib_bgzip_somatic {
// Runs bgzip on somatically tagged channels.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path(inf) - Input File 
//   path inf -  Input File
//
// output:
//   tuple => emit: bgzip_files
//     val(pat_name) - Patient Name
//     val(norm_run) - Normal Run Name
//     val(tumor_run) - Tumor Run Name
//     val(dataset) - Dataset
//     path("*gz") - bgzipped Output File
//
// require:
//   INPUT_FILE

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'htslib_container'
  label 'htslib_bgzip'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/htslib_bgzip"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(inf)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*gz"), emit: bgzip_files

  script:
  """
  bgzip ${inf}
  """
}


process htslib_bgzip {
// Runs bgzip
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(inf) - Input File 
//   path inf -  Input File
//
// output:
//   tuple => emit: bgzip_files
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*${inf}.gz") - bgzipped Output File
//
// require:
//   INPUT_FILE

  tag "${inf}"
  label 'htslib_container'
  label 'htslib_bgzip'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/htslib_bgzip"

  input:
  tuple val(pat_name), val(run), val(dataset), path(inf)

  output:
  tuple val(pat_name), val(run), val(dataset), path("${inf}.gz"), emit: bgzip_files

  script:
  """
  bgzip ${inf}
  """
}


process htslib_tabix {
// Runs tabix
//
// input:
//   path inf -  Input File
//
// output:
//   tuple => emit: tbx_file
//     path(inf) -  Input File
//     path("*tbi") - Tabix File
//
// require:
//   INPUT_FILE

  tag "${inf}"
  label 'htslib_container'
  label 'htslib_tabix'

  input:
  path inf

  output:
  tuple path(inf), path("*tbi"), emit: tabix_file

  script:
  def inf_proxy  = inf.name != 'dummy_file' ? "tabix -f ${inf}" : "touch ${inf}.tbi"
  """
  ${inf_proxy}
  """
}


process htslib_bgzip_ref {
// Runs bgzip
//
// input:
//   path inf -  Input File
//
// output:
//   tuple => emit: bgzip_files
//     path(inf) -  Input File
//     path("${inf}.gz") - Compressed Output File
//
// require:
//   INPUT_FILE

  tag "${inf}"
  label 'htslib_container'
  label 'htslib_bgzip_ref'

  input:
  path inf

  output:
  path("${inf}.gz"), emit: bgzip_files

  script:
  """
  bgzip ${inf}
  """
}
