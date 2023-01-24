#!/usr/bin/env nextflow

process sniffles {
// Runs sniffles
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

  tag "${dataset}/${pat_name}/${run}"
  label 'sniffles_container'
  label 'sniffles_mem'
  cache 'lenient'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam), path(bai)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path('*.vcf'), emit: vcfs

  script:
  """
  sniffles \
    --input ${bam} \
    --vcf ${dataset}-${pat_name}-${run}.sniffles.vcf \
    ${parstr}
  """
}
