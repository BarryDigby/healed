#!/usr/bin/env nextflow

include { samtools_view } from '../samtools/samtools.nf'

process hervquant_filter_alignments_sub {
// Subprocess responsible for filtering alignment files.

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*filt.sam"), emit: filt_sams

  script:
  """
  ALN_BUFFER="${aln}"
  ALN_FILT=\${ALN_BUFFER%.sam}.filt.sam

  sed '/uc.*/d' ${aln} > \${ALN_FILT}
  """
}


workflow hervquant_filter_alignments {
// Runs sed and samtools per HervQuant requirements.
// This should likely be broken up to be more consistent with other viral
// modules (virdetect).
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File
//
// output:
//   tuple => filt_sams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*filt.sam") - Filtered BAM File
//
// require:
//   BAMS
  take:
    bams
  main:
    hervquant_filter_alignments_sub(
      bams)
    samtools_view(
      hervquant_filter_alignments_sub.out.filt_bams,
      '-bS')
  emit:
    bams = samtools_view.out.bams
}
