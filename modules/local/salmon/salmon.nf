#!/usr/bin/env nextflow

process salmon_index {
// Runs salmon index command
//
// input:
//   path fa - Reference FASTA
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: idx_files
//     path(fa) - Reference FASTA
//     path("${fa}*index") - Index Files

// require:
//   params.salmon$reference
//   params.salmon$salmon_index_parameters

  tag "${fa}"
  label 'salmon_container'
  label 'salmon_index'
  storeDir "${params.shared_dir}/${fa}/salmon_index/"

  input:
  path fa
  val parstr

  output:
  tuple path("${fa}"), path("${fa}*index"), emit: idx_files

  script:
  """
  salmon index --threads ${task.cpus} ${parstr} -t ${fa} -i ${fa}.index
  """
}


process salmon_map_quant {
// Runs salmon quant (for quantifying from FASTQs)
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   tuple
//     path(fa) - Reference FA
//     path(idx_files) - Index Files
//   val parstr - Additoinal Parameters
//
// output:
//   tuple => emit: quants
//       val(pat_name) - Patient Name
//       val(run) - Run Name
//       val(dataset) - Dataset
//       path('*quant.sf') - Quant File

// require:
//   FQS
//   IDX_FILES
//   params.salmon$salmon_map_quant_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'salmon_container'
  label 'salmon_map'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/salmon_map_quant"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("${dataset}-${pat_name}-${run}.quant.sf"), emit: quants

  script:
  """
  salmon quant --threads ${task.cpus} -i ${idx} -l a -1 ${fq1} -2 ${fq2} -o . ${parstr}
  mv quant.sf ${dataset}-${pat_name}-${run}.quant.sf
  """
}



process salmon_aln_quant {
// Runs salmon quant (for quantifying from alignments)
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - FASTQ Prefix
//     val(dataset) - Dataset
//     path(aln) - Alignment File
//   path(fa) - Reference FA
//   val out_dir - Output Directory
//   val shared_dir - Shared Output Directory
//
// output:
//   tuple => emit: quants
//       val(pat_name) - Patient Name
//       val(run) - FASTQ Prefix
//       val(dataset) - Dataset
//       path('*quant.sf') - Quant File
//   path('salmon-*') - for publishDir

// require:
//   ALNS
//   REF
//   params.salmon$salmon_aln_quant_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'salmon_container'
  label 'salmon_aln'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/salmon_aln_quant"

  input:
  tuple val(pat_name), val(run), val(dataset), path(aln)
  path fa
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*sf'), emit: quants

  script:
  """
  salmon quant --threads ${task.cpus} -t ${fa}  -l a -a ${aln} -o . ${parstr}
  mv quant.sf ${dataset}-${pat_name}-${run}.quant.sf
  """
}
