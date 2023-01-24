#!/usr/bin/env nextflow

process awk_get_unaligned {
// required:
//   ALNS

  tag "${dataset}/${pat_name}/${run}"
  label 'virdetect_container'
//  publishDir "${params.samps_out_dir}/${dataset}/${run}/awk_get_unaligned"

  input:
  tuple val(pat_name), val(run), val(dataset), path(sam)

  output:
  tuple val(pat_name), val(run), val(dataset), path("unaligned.bam"), emit: bams

  script:
  """
  samtools view ${sam} | sh /opt/virdetect/awk_column3_star.sh | samtools view -bS - > unaligned.bam
  """
}


process awk_unaligned_bam_to_fastqs {
// Runs awk to extract reads from alignment file.
// required:
//   ALNS

  tag "${dataset}/${pat_name}/${run}"
  label 'virdetect_container'
//  publishDir "${params.samps_out_dir}/${dataset}/${run}/awk_unaligned_bam_to_fastqs"

  input:
  tuple val(pat_name), val(run), val(dataset), path(unaligned_bam)

  output:
  tuple val(pat_name), val(run), val(dataset), path("unaligned_1.fastq.gz"), path("unaligned_2.fastq.gz"), emit: fqs

  script:
  """
  samtools view ${unaligned_bam} | sh /opt/virdetect/awk_unalignedfq_1.sh > unaligned_1.fastq.tmp; samtools view ${unaligned_bam} | sh /opt/virdetect/awk_unalignedfq_2.sh > unaligned_2.fastq.tmp; python /opt/virdetect/only_paired_reads.py -1 unaligned_1.fastq.tmp -2 unaligned_2.fastq.tmp -o1 unaligned_1.fastq -o2 unaligned_2.fastq; gzip *fastq
  """
}



process count_star_viral_alignments {
// Runs countStarViralAlignments to determine which reads align to viral
// reference.

// require:
//   ALNS

  tag "${dataset}/${pat_name}/${run}"
  label 'virdetect_container'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/count_star_viral_alignments"

  input:
  tuple val(pat_name), val(run), val(dataset), path(sam)

  output:
  tuple val(pat_name), val(run), val(dataset), path('*viralReadCounts.txt'), emit: viral_counts

  script:
  """
  java -Xmx4G -cp /opt/virdetect/java/java-src/picard-1.92.jar:/opt/virdetect/java/java-src/sam-1.92.jar:/opt/virdetect/java/countStarViralAlignments.jar countStarViralAlignments ${dataset}-${pat_name}-${run} ${sam} ${dataset}-${pat_name}-${run}.viralReadCounts.txt
  """
}
