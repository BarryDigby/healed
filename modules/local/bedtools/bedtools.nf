#!/usr/bin/env nextflow

process bedtools_intersect {
// Compares two or more BED/BAM/VCF/GFF files and identifies all the regions in
// the genome where the features in the two files overlap (that is, share at
// least one base pair in common).
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File
//  path(bed) - BED File
//
// output:
//   tuple => emit: filt_bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*filt.bam")

// require:
//   ALNS
//   params.bedtools$bedtools_intersect_bed

  tag "${dataset}/${pat_name}/${run}"
  label 'bedtools_container'
  label 'bedtools_intersect'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  path bed

  output:
  tuple val(pat_name), val(run), val(dataset), path("*filt.bam"), emit: filt_bams

  script:
  """
  BAM_NAME=\$(echo ${bam})
  OUTPUT_NAME=\$(echo \${BAM_NAME%.bam})

  bedtools intersect -a ${bam} -b ${bed}  > \${OUTPUT_NAME}.filt.bam
  """
}


process bedtools_intersect_vcfs {
// Compares two or more BED/BAM/VCF/GFF files and identifies all the regions in
// the genome where the features in the two files overlap (that is, share at
// least one base pair in common).
//
// This is a specific use-case for LENS. This should ideally be wrapped up into
// the standard bedtools_intersect process.
//
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(bam) - BAM File
//  path(bed) - BED File
//
// output:
//   tuple => emit: filt_bams
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*filt.bam")

// require:
//   ALNS
//   params.bedtools$bedtools_intersect_bed

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/bedtools_intersect"
  label 'bedtools_container'
  label 'bedtools_intersect'

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcfs)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*intersect.vcf"), emit: intersect_vcfs

  script:
  """
  bedtools intersect -wa -u -header -a *mutect*vcf.gz -b *cadabra*vcf.gz -b *strelka*snv*vcf.gz -b *strelka*indel*vcf.gz > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.intersect.vcf
  """
}
