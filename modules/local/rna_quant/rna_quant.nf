#!/usr/bin/env nextflow

include { abra2_rna } from '../abra2/abra2.nf'
include { bedtools_intersect } from '../bedtools/bedtools.nf'
include { gffread_make_tx_fa } from '../gffread/gffread.nf'
include { kallisto_index }     from '../kallisto/kallisto.nf'
include { kallisto_quant }     from '../kallisto/kallisto.nf'
include { salmon_index }       from '../salmon/salmon.nf'
include { salmon_aln_quant }   from '../salmon/salmon.nf'
include { salmon_map_quant }   from '../salmon/salmon.nf'
include { star_index }         from '../star/star.nf'
include { star_map }           from '../star/star.nf'
include { trim_galore }        from '../trim_galore/trim_galore.nf'
include { get_fastqs }         from '../utilities/utilities.nf'
include { samtools_faidx }     from '../samtools/samtools.nf'
include { samtools_index }     from '../samtools/samtools.nf'
include { picard_create_seq_dict } from '../picard2/picard2.nf'

workflow manifest_to_salmon_counts {
// require:
//   params.rna_quant$manifest_to_salmon_counts$salmon_ref
//   MANIFEST
  take:
    salmon_ref
    manifest
  main:
    get_fastqs(
      manifest,
      params.fq_dir)
    raw_fqs_to_salmon_counts(
      salmon_ref,
      get_fastqs.out.fastqs)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_salmon_counts.out.procd_fqs
    quants = raw_fqs_to_salmon_counts.out.quants
}


workflow manifest_to_kallisto_counts {
// require:
//   params.rna_quant$manifest_to_kallisto_counts$kallisto_ref
//   MANIFEST
  take:
    kallisto_ref
    manifest
  main:
    get_fastqs(
      manifest,
      params.rna_quant$fq_dir)
    raw_fqs_to_kallisto_counts(
      kallisto_ref,
      get_fastqs.out.fastqs)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_kallisto_counts.out.procd_fqs
    quants = raw_fqs_to_kallisto_counts.out.quants
}


workflow manifest_to_star_alns {
// require:
//   params.rna_quant$manifest_to_star_alns$star_ref
//   MANIFEST
  take:
    star_ref
    manifest
  main:
    get_fastqs(
      manifest,
      params.rna_quant$fq_dir)
    raw_fqs_to_star_alns(
      star_ref,
      get_fastqs.out.fastqs)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_star_alns.out.procd_fqs
    alns = raw_fqs_to_star_alns.out.alns
    alt_alns = raw_fqs_to_star_alns.out.alt_alns
}


workflow manifest_to_indel_realigned_star_alns {
// require:
//   params.rna_quant$manifest_to_indel_realigned_star_alns$star_ref
//   params.rna_quant$manifest_to_indel_realigned_star_alns$gtf
//   MANIFEST
  take:
    star_ref
    gtf
    manifest
  main:
    get_fastqs(
      manifest,
      params.rna_quant$fq_dir)
    raw_fqs_to_indel_realigned_star_alns(
      star_ref,
      gtf,
      get_fastqs.out.fastqs)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_star_alns.out.procd_fqs
    alns = raw_fqs_to_star_alns.out.alns
    alt_alns = raw_fqs_to_star_alns.out.alt_alns
}


workflow manifest_to_star_alns_salmon_counts {
// require:
//   params.rna_quant$manifest_to_star_alns_salmon_counts$star_salmon_ref
//   params.rna_quant$manifest_to_star_alns_salmon_counts$gtf
//   MANIFEST
//   params.rna_quant$manifest_to_star_alns_salmon_counts$trim_galore_parameters
//   params.rna_quant$manifest_to_star_alns_salmon_counts$star_index_parameters
//   params.rna_quant$manifest_to_star_alns_salmon_counts$star_map_parameters
//   params.rna_quant$manifest_to_star_alns_salmon_counts$star_alt_capture
//   params.rna_quant$manifest_to_star_alns_salmon_counts$salmon_aln_quant_parameters
  take:
    star_salmon_ref
    gtf
    manifest
    trim_galore_parameters
    star_index_parameters
    star_map_parameters
    star_alt_capture
    salmon_aln_quant_parameters 
  main:
    get_fastqs(
      manifest,
      params.fq_dir)
    raw_fqs_to_star_alns_salmon_counts(
      star_salmon_ref,
      gtf,
      get_fastqs.out.fastqs,
      trim_galore_parameters,
      star_index_parameters,
      star_map_parameters,
      star_alt_capture,
      salmon_aln_quant_parameters)
  emit:
    fqs = get_fastqs.out.fastqs
    procd_fqs = raw_fqs_to_star_alns_salmon_counts.out.procd_fqs
    alns = raw_fqs_to_star_alns_salmon_counts.out.alns
    alt_alns = raw_fqs_to_star_alns_salmon_counts.out.alt_alns
    junctions = raw_fqs_to_star_alns_salmon_counts.out.junctions
    quants = raw_fqs_to_star_alns_salmon_counts.out.quants
    star_logs = raw_fqs_to_star_alns_salmon_counts.out.star_logs
}


workflow raw_fqs_to_salmon_counts {
// require:
//   params.rna_quant$raw_fqs_to_salmon_counts$salmon_ref
//   FQS
  take:
    salmon_ref
    fqs
  main:
    trim_galore(
      fqs,
      params.rna_quant$raw_fqs_to_salmon_counts$trim_galore_parameters)
    procd_fqs_to_salmon_counts(
      salmon_ref,
      trim_galore.out.procd_fqs)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    quants = procd_fqs_to_salmon_counts.out.quants
}


workflow raw_fqs_to_kallisto_counts {
// require:
//   params.rna_quant$raw_fqs_to_kallisto_counts$kallisto_ref
//   params.rna_quant$raw_fqs_to_kallisto_counts$fqs
  take:
    kallisto_ref
    fqs
  main:
    trim_galore(
      fqs,
      params.rna_quant$raw_fqs_to_kallisto_counts$trim_galore_parameters,
      params.rna_quant$samps_out_dir)
    procd_fqs_to_kallisto_counts(
      kallisto_ref,
      trim_galore.out.procd_fqs)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    quants = procd_fqs_to_kallisto_counts.out.quants
}


workflow raw_fqs_to_star_alns {
// require:
//   params.rna_quant$raw_fqs_to_star_alns$star_ref
//   FQS
  take:
    star_ref
    fqs
  main:
    trim_galore(
      fqs,
      params.rna_quant$raw_fqs_to_star_alns$trim_galore_parameters)
    procd_fqs_to_star_alns(
      star_ref,
      trim_galore.out.procd_fqs)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    alns = procd_fqs_to_star_alns.out.alns
    alt_alns = procd_fqs_to_star_alns.out.alt_alns
    junctions = procd_fqs_to_star_alns.out.junctions
}


workflow raw_fqs_to_indel_realigned_star_alns {
// require:
//   params.rna_quant$raw_fqs_to_indel_realigned_star_alns$star_ref
//   params.rna_quant$raw_fqs_to_indel_realigned_star_alns$gtf
//   FQS
  take:
    star_ref
    gtf
    fqs
  main:
    trim_galore(
      fqs,
      params.rna_quant$raw_fqs_to_star_alns$trim_galore_parameters)
    procd_fqs_to_indel_realigned_star_alns(
      star_ref,
      gtf,
      trim_galore.out.procd_fqs)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    alns = procd_fqs_to_star_alns.out.alns
    alt_alns = procd_fqs_to_star_alns.out.alt_alns
    juncitons = procd_fqs_to_star_alns.out.junctions
}


workflow raw_fqs_to_star_alns_salmon_counts {
// require:
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$star_salmon_ref
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$gtf
//   FQS
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$trim_galore_parameters
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$star_index_parameters
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$star_map_parameters
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$star_alt_capture
//   params.rna_quant$raw_fqs_to_star_alns_salmon_counts$salmon_aln_quant_parameters
  take:
    star_salmon_ref
    gtf
    fqs
    trim_galore_parameters
    star_index_parameters
    star_map_parameters
    star_alt_capture
    salmon_aln_quant_parameters
  main:
    trim_galore(
      fqs,
      trim_galore_parameters)
    procd_fqs_to_star_alns_salmon_counts(
      star_salmon_ref,
      gtf,
      trim_galore.out.procd_fqs,
      star_index_parameters,
      star_map_parameters,
      star_alt_capture,
      salmon_aln_quant_parameters)
  emit:
    procd_fqs = trim_galore.out.procd_fqs
    alns = procd_fqs_to_star_alns_salmon_counts.out.alns
    alt_alns = procd_fqs_to_star_alns_salmon_counts.out.alt_alns
    junctions = procd_fqs_to_star_alns_salmon_counts.out.junctions
    quants = procd_fqs_to_star_alns_salmon_counts.out.quants
    star_logs = procd_fqs_to_star_alns_salmon_counts.out.star_logs
}


workflow procd_fqs_to_salmon_counts {
// require:
//   params.rna_quant$procd_fqs_to_salmon_counts$salmon_ref
//   FQS
  take:
    salmon_ref
    fqs
  main:
    salmon_index(
      salmon_ref,
      params.rna_quant$procd_fqs_to_salmon_counts$salmon_index_parameters)
    salmon_map_quant(
      fqs,
      salmon_index.out.idx_files,
      params.rna_quant$procd_fqs_to_salmon_counts$salmon_map_quant_parameters)
  emit:
    quants = salmon_map_quant.out.quants
}


workflow procd_fqs_to_kallisto_counts {
// require:
//   params.rna_quant$procd_fqs_to_kallisto_counts$kallisto_ref
//   params.rna_quant$procd_fqs_to_kallisto_counts$fqs
  take:
    kallisto_ref
    fqs
  main:
    kallisto_index(
      kallisto_ref,
      params.rna_quant$procd_fqs_to_kallisto_counts$kallisto_index_parameters,
      params.rna_quant$indices_out_dir,
      params.shared_dir)
    kallisto_quant(
      fqs,
      kallisto_index.out.idx_files,
      params.rna_quant$procd_fqs_to_kallisto_counts$kallisto_quant_parameters,
      params.rna_quant$samps_out_dir,
      params.shared_dir)
  emit:
    quants = kallisto_quant.out.quants
}


workflow procd_fqs_to_star_alns {
// require:
//   params.rna_quant$procd_fqs_to_star_alns$star_ref
//   FQS
  take:
    star_ref
    fqs
    gtf
    star_index_parameters
    star_map_parameters
    star_alt_capture
  main:
    star_index(
      star_ref,
      star_index_parameters)
    star_map(
      fqs,
      star_index.out.idx_files,
      star_map_parameters,
      star_alt_capture,
      gtf)
  emit:
    alns = star_map.out.alns
    alt_alns = star_map.out.alt_alns
    junctions = star_map.out.junctions
    standard_junctions = star_map.out.standard_junctions
    star_logs = star_map.out.star_logs
}


workflow procd_fqs_to_indel_realigned_star_alns {
// require:
//   params.rna_quant$procd_fqs_to_indel_realigned_star_alns$star_ref
//   params.rna_quant$procd_fqs_to_indel_realigned_star_alns$gtf
//   FQS
  take:
    star_ref
    gtf
    fqs
    star_index_parameters
    star_map_parameters
    star_alt_capture
    samtools_index_parameters
    abra2_parameters
  main:
    star_index(
      star_ref,
      star_index_parameters)
    star_map(
      fqs,
      star_index.out.idx_files,
      star_map_parameters,
      star_alt_capture,
      gtf)
    bedtools_intersect(
      star_map.out.alns,
      targets_bed)
    make_ancillary_index_files(
      star_ref)
    samtools_index(
      bedtools_intersect.out.filt_alns,
      samtools_index_parameters)
    abra2_rna(
      samtools_index.out.bams_and_bais,
      make_ancillary_index_files.out.collective_idx_files,
      gtf,
      abra2_parameters,
      'tmp_dir')
  emit:
    abra_alns = abra2_rna.out.abra_bams
    alns = star_map.out.alns
    alt_alns = star_map.out.alt_alns
}


workflow procd_fqs_to_star_alns_salmon_counts {
// require:
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$star_salmon_ref
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$gtf
//   FQS
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$star_index_parameters
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$star_map_parameters
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$star_alt_capture
//   params.rna_quant$procd_fqs_to_star_alns_salmon_counts$salmon_aln_quant_parameters
  take:
    star_salmon_ref
    gtf
    fqs
    star_index_parameters
    star_map_parameters
    star_alt_capture
    salmon_aln_quant_parameters
  main:
    salmon_ref = star_salmon_ref
    procd_fqs_to_star_alns(
      star_salmon_ref,
      fqs,
      gtf,
      star_index_parameters,
      star_map_parameters,
      star_alt_capture)
    //If gtf, then create txome fa and set it to salmon_ref
    if( gtf )
      gffread_make_tx_fa(
        star_salmon_ref,
        gtf)
      gffread_make_tx_fa.out.tx_fa.set{ salmon_ref }
    //If alt_alns exists, then use it. Otherwise, defer to standard alns.
    if( procd_fqs_to_star_alns.out.alt_alns )
      star_alns_to_salmon_counts(
        procd_fqs_to_star_alns.out.alt_alns,
        salmon_ref,
        salmon_aln_quant_parameters)
    else
      star_alns_to_salmon_counts(
        procd_fqs_to_star_alns.out.alns,
        salmon_ref,
        salmon_aln_quant_parameters)
  emit:
    alns = procd_fqs_to_star_alns.out.alns
    alt_alns = procd_fqs_to_star_alns.out.alt_alns
    junctions = procd_fqs_to_star_alns.out.junctions
    standard_junctions = procd_fqs_to_star_alns.out.standard_junctions
    quants = star_alns_to_salmon_counts.out.quants
    star_logs = procd_fqs_to_star_alns.out.star_logs
}


workflow star_alns_to_salmon_counts {
  take:
    alns
    salmon_ref
    salmon_aln_quant_parameters
  main:
    salmon_aln_quant(
      alns,
      salmon_ref,
      salmon_aln_quant_parameters)
  emit:
    quants = salmon_aln_quant.out.quants
}


workflow make_ancillary_index_files {
  take:
    ref
  main:
    samtools_faidx(
      ref,
      '')
//      params.make_ancillary_index_files$samtools_faidx_parameters)
    picard_create_seq_dict(
      ref,
      '')
    //Can't join since the ref symlinks are different.
    samtools_faidx.out.faidx_file
      .concat(picard_create_seq_dict.out.dict_file)
      .collect().map{ [it[0], it[1], it[3]] }
      .set{ collective_idx_files }
  emit:
    collective_idx_files
}
