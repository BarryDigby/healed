#!/usr/bin/env nextflow

include { bwa_index } from '../bwa/bwa.nf'
include { bwa_mem_samtools_sort } from '../bwa/bwa.nf'

include { manifest_to_dna_procd_fqs } from '../preproc/preproc.nf'

include { samtools_index } from '../samtools/samtools.nf'
include { samtools_sort } from '../samtools/samtools.nf'
include { samtools_view } from '../samtools/samtools.nf'

include { sequenza_gc_wiggle } from '../sequenza/sequenza.nf'
include { sequenza_bam2seqz } from '../sequenza/sequenza.nf'
include { sequenza_seqz_binning } from '../sequenza/sequenza.nf'
include { sequenza_extract } from '../sequenza/sequenza.nf'
include { sequenza_fit }  from '../sequenza/sequenza.nf'
include { sequenza_result } from '../sequenza/sequenza.nf'

include { starfusion } from '../starfusion/starfusion.nf'

workflow manifest_to_sequenza {
// require:
//   params.onco$manifest_to_sequenza$dna_ref
//   params.onco$manifest_to_sequenza$manifest
  take:
    dna_ref
    manifest
  main:
    manifest_to_dna_procd_fqs(
      manifest)
    procd_fqs_to_sequenza(
      dna_ref,
      manifest,
      manifest_to_dna_procd_fqs.out.procd_fqs)
}

//workflow raw_fqs_to_sequenza {
//  take:
//  main:
//  emit:
//}

//workflow raw_alnz_to_sequenza {
//  take:
//    ref
//    manifest
//    bams
//    bais
//  main:
//    ref.view()
//    bams.view()
//    bais.view()
//}

workflow procd_fqs_to_sequenza {
  take:
    dna_ref
    manifest
    fqs
  main:
    bwa_index(
      dna_ref,
      params.onco$procd_fqs_to_sequenza$bwa_index_params)
    bwa_mem_samtools_sort(
      fqs,
      bwa_index.out.idx_files,
      params.onco$procd_fqs_sequenza$bwa_mem_params)
//    samtools_view(
//      bwa_mem.out.sams,
//      params.onco$procd_fqs_to_sequenza$samtools_view)
//    samtools_sort(
//      bwa_mem_samtools_view.out.bams,
//      '')
    samtools_index(
      bwa_mem_samtools_sort.out.bams,
      '')
//    raw_alnz_to_sequenza(dna_ref,
//                         manifest,
//                         samtools_sort.out.bams,
//                         samtools_index.out.bais)
    //A lot of filtering, joining, etc. to ensure inputs go in as expected.
    manifest.filter{ it[5] =~ /TRUE/ }.set{ norms }
    manifest.filter{ it[5] =~ /FALSE/ }.set{ tumors }
    samtools_sort.out.bams.join(samtools_index.out.bais, by: [0, 1]).set{ bams_bais }
    bams_bais.join(norms, by: [0, 1]).map{ [it[0], it[1], it[2], it[3], it[5]] }.set{ norm_bams_bais }
    bams_bais.join(tumors, by: [0, 1]).map{ [it[0], it[1], it[2], it[3], it[5]] }.set{ tumor_bams_bais }
    norm_bams_bais.join(tumor_bams_bais, by: [0, 2]).set{ norm_tumor_bams_bais }
    sequenza_gc_wiggle(dna_ref, '50', '')
    sequenza_bam2seqz(norm_tumor_bams_bais, sequenza_gc_wiggle.out.gc_wig)
    sequenza_seqz_binning(sequenza_bam2seqz.out.seqz, '50')
    sequenza_extract(sequenza_seqz_binning.out.small_seqz)
    sequenza_fit(sequenza_extract.out.extracts)
//    sequenza_extract.out.extracts.join(sequenza_fit.out.fits, by: 0).set{ collective_files }
    sequenza_result(sequenza_fit.out.fits)
}


workflow procd_alns_to_sequenza {
  take:
    paired_alns
    dna_ref
    bed
    sequenza_gc_wiggle_window_size
    sequenza_gc_wiggle_parameters
    sequenza_seqz_binning_window_size
  main:
    //A lot of filtering, joining, etc. to ensure inputs go in as expected.
    sequenza_gc_wiggle(
      dna_ref,
      sequenza_gc_wiggle_window_size,
      sequenza_gc_wiggle_parameters)
    sequenza_bam2seqz(
       paired_alns,
       sequenza_gc_wiggle.out.gc_wig,
       bed)
    sequenza_seqz_binning(
      sequenza_bam2seqz.out.seqz,
      sequenza_seqz_binning_window_size)
    sequenza_extract(
      sequenza_seqz_binning.out.small_seqz)
    sequenza_fit(
      sequenza_extract.out.extracts)
    sequenza_result(sequenza_fit.out.fits)
  emit:
    segments = sequenza_result.out.segments 
    segments_and_solutions = sequenza_result.out.segments_and_solutions
}
