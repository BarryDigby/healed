#!/usr/bin/env nextflow

process gffread_make_tx_fa {
// GFF/GTF utility providing format conversions, filtering, FASTA sequence
// extraction and more.
//
// input:
//   path fa - Reference FASTA
//   path gtf - Reference GTF
//
// output:
//   path("transcripts.fa") => emit: tx_fa

// require:
//   params.gffread$gffread_make_tx_fa$dna_ref
//   params.gffread$gffread_make_tx_fa$gtf

  tag "${fa}/${gtf}"
  label 'gffread_container'
  label 'gffread'

  input:
  path fa
  path gtf

  output:
  path "*transcripts.fa", emit: tx_fa

  script:
  """
  AFA=\$(echo ${fa})
  AGTF=\$(echo ${gtf})
  gffread -w \${AFA%.fa*}.\${AGTF%.gtf}.transcripts.fa -g ${fa} ${gtf}
  """
}
