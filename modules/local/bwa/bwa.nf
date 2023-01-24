#!/usr/bin/env nextflow

process bwa_index {
// Runs bwa index
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

  storeDir "${params.shared_dir}/${fa}/bwa_index"
  tag "${fa}"
  label 'bwa_container'
  label 'bwa_index'

  input:
  path fa
  val parstr

  output:
  tuple path(fa), path("${fa}.*"), emit: idx_files

  script:
  """
  bwa index ${parstr} ${fa}
  """
}


process bwa_mem {
// Runs bwa mem
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
  label 'bwa_container'
  label 'bwa_mem'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx_files)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.sam'), emit: sams

  script:
  """
  bwa mem ${parstr} ${fa} ${fq1} ${fq2} > ${dataset}-${pat_name}-${run}.sam -t ${task.cpus}
  """
}


process bwa_mem_samtools_sort {
// Runs bwa mem piped to samtools sort (to minimize storage footprint)
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
//   tuple => emit: bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*.bam") - Output BAM File

// require:
//   FQS
//   IDX_FILES
//   params.bwa$bwa_mem_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'bwa_samtools_container'
  label 'bwa_samtools'

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx_files)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*sorted.bam'), emit: bams

  script:
  """
  bwa mem ${parstr} -R "@RG\\tID:${dataset}-${pat_name}-${run}\\tSM:${dataset}-${pat_name}-${run}\\tLB:NULL\\tPL:Illumina" ${fa} ${fq1} ${fq2} -t ${task.cpus} > tmp.bam
  samtools sort -@ ${task.cpus} tmp.bam > ${dataset}-${pat_name}-${run}.sorted.bam
  rm tmp.bam
  """
}
