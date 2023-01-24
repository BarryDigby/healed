#!/usr/bin/env nextflow

process snpeff_download {
// require:
//   REF

  tag "${ref}"
  label 'snpeff_container'
  label 'snpeff_download'

  input:
  val ref

  output:
  path("${ref}"), emit: snpeff_database

  script:
  """
  /opt/snpeff/snpeff/bin/snpEff download ${ref} -dataDir ${PWD}
  """
}



process snpeff_ann {
// require:
//   VCFS
//   params.snpeff$snpeff_ann_ref

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'snpeff_container'
  label 'snpeff_ann'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/snpeff_ann"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  path ref

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*annot.vcf"), emit: annot_vcfs

  script:
  """
  ORIG="${vcf}"
  /opt/snpeff/snpeff/bin/snpEff ann ${ref} ${vcf} -noShiftHgvs -c ${ref}/snpEff.config > \${ORIG%.vcf*}.annot.vcf -dataDir \${PWD}
  """
}


process snpeff_ann_germline {
// require:
//   VCFS
//   params.snpeff$snpeff_ann_germline_ref

  tag "${dataset}/${pat_name}/${run}"
  label 'snpeff_container'
  label 'snpeff_ann_germline'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/snpeff_ann"

  input:
  tuple val(pat_name), val(run), val(dataset), path(vcf)
  path ref

  output:
  tuple val(pat_name), val(run), val(dataset), path("*annot.vcf"), emit: annot_vcfs

  script:
  """
  ORIG="${vcf}"
  /opt/snpeff/snpeff/bin/snpEff ann ${ref} -dataDir \${PWD} ${vcf} > \${ORIG%.vcf*}.annot.vcf
  """
}


process snpsift_filter {
// require:
//   VCFS
//   params.snpeff$snpsift_filter_parameters

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'snpeff_container'
  label 'snpeff_filter'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/snpeff_filter"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf)
  val parstr
  val suffix

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*${suffix}.vcf"), optional: true, emit: filt_vcfs

  script:
  """
  ORIG="${vcf}"
  /opt/snpeff/snpeff/bin/SnpSift filter "${parstr}" -f ${vcf} > \${ORIG%.vcf*}.${suffix}.vcf
  if [ `grep -v "^#" \${ORIG%.vcf*}.${suffix}.vcf | wc -l` == 0 ]; then
    mv \${ORIG%.vcf*}.${suffix}.vcf  \${ORIG%.vcf*}.${suffix}_no_hits.vcf
  fi
  """
}
