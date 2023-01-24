#!/usr/bin/env nextflow

process optitype_razers3 {
//Technically, this is its own separate aligner, but it's a highly suggested
//step for filtering to HLA reads.

  tag "${dataset}/${pat_name}/${prefix}"
  label 'optitype_container'
  label 'optitype_razers3'

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*_1.bam"), path("*_2.bam"), emit: filtd_bams

  script:
  """
  /usr/local/bin/razers3 -i 95 -m 1 -dr 0 -o ${dataset}-${pat_name}-${prefix}_1.bam /usr/local/bin/OptiType-1.3.4/data/hla_reference_dna.fasta ${fq1}&
  /usr/local/bin/razers3 -i 95 -m 1 -dr 0 -o ${dataset}-${pat_name}-${prefix}_2.bam /usr/local/bin/OptiType-1.3.4/data/hla_reference_dna.fasta ${fq2}&
  wait
  """
}


process optitype_samtools_bam2fq {

  tag "${dataset}/${pat_name}/${prefix}"

  label 'optitype_container'
  label 'optitype_samtools_bam2fq'

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(bam1), path(bam2)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*_1.fastq"), path("*_2.fastq"), emit: filtd_fqs

  script:
  """
  samtools bam2fq ${bam1} > ${dataset}-${pat_name}-${prefix}_1.fastq
  samtools bam2fq ${bam2} > ${dataset}-${pat_name}-${prefix}_2.fastq
  """
}

process optitype_sub {

  tag "${dataset}/${pat_name}/${prefix}"

  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${prefix}/optitype"
  label 'optitype_container_alt'
  label 'optitype'
  containerOptions '--cleanenv --no-home -B /datastore'

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*optitype_calls*"), emit: calls

  script:
  """
  python /usr/local/bin/OptiType/OptiTypePipeline.py -i ${fq1} ${fq2} --dna --outdir calls
  mv calls/*/*result.tsv ${dataset}-${pat_name}-${prefix}.optitype_calls.tsv
  """
}

workflow optitype {
  take:
    fqs
  main:
    optitype_razers3(
      fqs,
      '')
    optitype_samtools_bam2fq(
      optitype_razers3.out.filtd_bams,
      '')
    optitype_sub(
      optitype_samtools_bam2fq.out.filtd_fqs,
      '')
  emit:
   calls = optitype_sub.out.calls
}
