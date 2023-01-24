#!/usr/bin/env nextflow

process bbmap_index {
// Runs bwa index
//
// input:
//   path fa - Reference FASTA
//   val params - Additional Parameters
//   val out_dir - Output Directory
//   val shared_dir - Shared Output Directory
//
// output:
//   tuple => emit: idx_files
//     path(fa) - Reference FASTA
//     path("${fa}.*") - Index Files

// require:
//   params.bwa$dna_reference
//   params.bwa$bbmap_index_parameters

  label 'bbmap_container'
  label 'bbmap_index'
  tag "${fa}"
  storeDir "${params.shared_dir}/${fa}/bbmap_index"

  input:
  path fa
  val parstr

  output:
  tuple path(fa), path("ref/"), emit: idx_files

  script:
  """
  bbmap.sh -Xmx${task.memory.toGiga()}g ref=${fa} ${parstr}
  """
}


process bbmap {
// Runs bbmap
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
//   FQS
//   IDX_FILES
//   params.bbmap$bbmap_parameters

  label 'bbmap_container'
  label 'bbmap'
  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx_files)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.sam'), emit: alns

  script:
  """
  bbmap.sh -Xmx${task.memory.toGiga()}g ref=${fa} in=${fq1} in2=${fq2} out=${dataset}-${pat_name}-${run}.sam ${parstr}
  """
}

process bbmap_samtools_sort {
// Runs bbmap
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
//     path('*.bam') - Alignment Output File

// require:
//   FQS
//   IDX_FILES
//   params.bbmap$bbmap_parameters

  label 'bbmap_container'
  label 'bbmap_samtools'
  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx_files)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.bam'), emit: alns

  script:
  """
  bbmap.sh -Xmx${task.memory.toGiga()}g ref=${fa} in=${fq1} in2=${fq2} out=tmp.sam ${parstr} t=${task.cpus}
  samtools sort -@ ${task.cpus} -o bam tmp.sam  >  ${dataset}-${pat_name}-${run}.sorted.bam
  rm -rf tmp.sam
  """
}
