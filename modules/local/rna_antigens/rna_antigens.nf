#!/usr/bin/env nextflow

include { neosplice } from '../neosplice/neosplice.nf'


workflow bams_to_neosplice {
// require:
//   MANIFEST
//   RNA_BAMS
//   params.rna_antigens$bams_to_neosplice$neosplice_ref
//   params.rna_antigens$bams_to_neosplice$gff
//   HLAS
//   params.rna_antigens$bams_to_neosplice$neosplice_augmented_splice_graph_build_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_get_max_kmer_length_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_convert_bams_to_fasta_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_get_splice_junctions_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_msbwt_is_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_convert_bwt_format_parameters
//   params.rna_antigens$bams_to_neosplice$neosplcie_kmer_search_bwt_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_search_bam_parameters
//   params.rna_antigens$bams_to_neosplice$samtools_sort_parameters
//   params.rna_antigens$bams_to_neosplice$samtools_index_parameters
//   params.rna_antigens$bams_to_neosplice$neosplice_kmer_graph_inference_parameters
  take:
    manifest
    rna_bams
    fa
    gff
    hlas
    neosplice_augmented_splice_graph_build_parameters
    neosplice_get_max_kmer_length_parameters
    neosplice_convert_bams_to_fasta_parameters
    neosplice_get_splice_junctions_parameters
    neosplice_msbwt_is_parameters
    neosplice_convert_bwt_format_parameters
    neosplice_kmer_search_bwt_parameters
    neosplice_search_bam_parameters
    samtools_sort_parameters
    samtools_index_parametrs
    neosplice_kmer_graph_inference_parameters
  main:
    neosplice(
      manifest,
      rna_bams,
      fa,
      gff,
      hlas,
      neosplice_augmented_splice_graph_build_parameters,
      neosplice_get_max_kmer_length_parameters,
      neosplice_convert_bams_to_fasta_parameters,
      neosplice_get_splice_junctions_parameters,
      neosplice_msbwt_is_parameters,
      neosplice_convert_bwt_format_parameters,
      neosplice_kmer_search_bwt_parameters,
      neosplice_search_bam_parameters,
      samtools_sort_parameters,
      samtools_index_parametrs,
      neosplice_kmer_graph_inference_parameters
       )
  emit:
    neoantigen_results = neosplice.out.neoantigen_results
    //outcome_peptides = neosplice.out.outcome_peptides
}
