#!/usr/bin/env nextflow

process seqtk_sample {
// require:
//   FQS
//   params.seqtk$seqtk_sample_read_count
//   SEED
//   params.seqtk$seqtk_sample_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'seqtk_container'
  label 'seqtk_sample'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/seqtk_sample"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val count
  val seed
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*_1.*.subd.*.fastq.gz"), path("*_2.*.subd.*.fastq.gz"), emit: subd_fqs

  script:
  """
  FQ1_BUFFER=`echo ${fq1}`
  FQ2_BUFFER=`echo ${fq2}`

  FQ1_OUT=`echo \${FQ1_BUFFER} | sed 's/.fastq.gz//g' | sed 's/.fq.gz//g'`
  FQ2_OUT=`echo \${FQ2_BUFFER} | sed 's/.fastq.gz//g' | sed 's/.fq.gz//g'`

  HR=`numfmt --to=si ${count}`

  seqtk sample -s${seed} ${parstr} ${fq1} ${count} | gzip -c > \${FQ1_OUT}${suffix}.subd.\${HR}.fastq.gz&
  seqtk sample -s${seed} ${parstr} ${fq2} ${count} | gzip -c > \${FQ2_OUT}${suffix}.subd.\${HR}.fastq.gz
  """
}

process seqtk_subseq {
// require:
//   FQS_AND_READ_NAMES
//   params.seqtk$seqtk_subseq_parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'seqtk_container'
  label 'seqtk_subseq'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/seqtk_subseq"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2), path(read_names)
  val suffix
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*extracted.fastq"), emit: extracted_fqs

  script:
  """
  seqtk subseq <(zcat ${fq1}) ${read_names}  > ${dataset}-${pat_name}-${run}${suffix}.extracted.fastq
  seqtk subseq <(zcat ${fq2}) ${read_names}  >> ${dataset}-${pat_name}-${run}${suffix}.extracted.fastq
  """
}
