#!/usr/bin/env nextflow

process minimap2_index {
// Creates a minimap2 index
//
// input:
//   path fa - Reference FASTA
//   val params - Additional Parameters
//
// output:
//   tuple => emit: idx_files
//     path(fa) - Reference FASTA
//     path("${fa}.*") - Index Files

// require:
//   params.bwa$dna_ref
//   params.bwa$bwa_index_parameters

  storeDir "${params.shared_dir}/${fa}/minimap2_index"
  tag "${fa}"
  label 'minimap2_container'
  label 'minimap2_index'

  input:
  path fa
  val parstr

  output:
  tuple path(fa), path("${fa}.*"), emit: idx_files

  script:
  """
#  minimap2 index ${parstr} ${fa}
  """
}


process minimap2 {
// Runs minimap2
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
//     path('*.sam') - Output SAM File

// require:
//   FQS
//   IDX_FILES
//   params.bwa$bwa_mem_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'minimap2_container'
  label 'minimap2'
  cache 'lenient'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq)
  tuple path(fa)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.paf'), optional: true, emit: pafs
  tuple val(pat_name), val(run), val(dataset), path('*.sam'), optional: true, emit: sams

  script:
  """
  minimap2 ${parstr} ${fa} ${fq} -t ${task.cpus} > ${dataset}-${pat_name}-${run}.sam
  """
}


process minimap2_samtools_sort {
// Runs minimap2
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
//     path('*.sam') - Output SAM File

// require:
//   FQS
//   IDX_FILES
//   params.bwa$bwa_mem_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'minimap2_samtools_container'
  label 'minimap2'
  cache 'lenient'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq)
  tuple path(fa)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.sorted.bam'), optional: true, emit: bams

  script:
  """
  minimap2 ${parstr} ${fa} ${fq} -t ${task.cpus - 1} |\
  samtools sort -@ 1 -o ${dataset}-${pat_name}-${run}.sorted.bam -
  """
}
