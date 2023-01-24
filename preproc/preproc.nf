#!/usr/bin/env nextflow

include { trim_galore } from '../trim_galore/trim_galore.nf'
include { get_fastqs } from '../utilities/utilities.nf'

workflow manifest_to_raw_fqs {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest.map{ [it[0], it[1], it[2], it[3]] }, //tuple val(pat_name), val(prefix), val(dataset), val(run)
      params.fq_dir)
  emit:
    fastqs = get_fastqs.out.fastqs
}


workflow manifest_to_rna_procd_fqs {
// require:
//   MANIFEST
//   params.manifest_to_rna_procd_fqs$trim_galore_parameters
  take:
    manifest
    trim_galore_parameters
  main:
    manifest_to_raw_fqs(
      manifest.filter{ it[4] =~ /RNA/ })
    trim_galore(
      manifest_to_raw_fqs.out.fastqs,
      trim_galore_parameters)
  emit:
    raw_fqs = manifest_to_raw_fqs.out.fastqs
    procd_fqs = trim_galore.out.procd_fqs
    fastqc_zips = trim_galore.out.fastqc_zips
}


workflow manifest_to_dna_procd_fqs {
// require:
//   MANIFEST
//   params.manifest_to_rna_procd_fqs$trim_galore_parameters
  take:
    manifest
    trim_galore_parameters
  main:
    manifest_to_raw_fqs(
      manifest.filter{ it[4] =~ /WES|WXS|DNA/ })
    trim_galore(
      manifest_to_raw_fqs.out.fastqs,
      trim_galore_parameters)
  emit:
    raw_fqs = manifest_to_raw_fqs.out.fastqs
    procd_fqs = trim_galore.out.procd_fqs
    fastqc_zips = trim_galore.out.fastqc_zips
}


workflow manifest_to_procd_fqs {
// require:
//   params.preproc$manifest_to_faw_fqs$manifest
  take:
    manifest
  main:
    manifest_to_raw_fqs(
      manifest)
    trim_galore(
      manifest_to_raw_fqs.out.fastqs,
      params.preproc$manifest_to_procd_fqs$trim_galore_parameters)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
}
