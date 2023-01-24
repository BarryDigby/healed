#!/usr/bin/env nextflow

include { lenstools_make_snv_peptides } from '../lenstools/lenstools.nf'
include { lenstools_make_snv_peptides_context } from '../lenstools/lenstools.nf'
include { lenstools_filter_mutant_snv_peptides } from '../lenstools/lenstools.nf'
include { lenstools_filter_mutant_indel_peptides } from '../lenstools/lenstools.nf'
include { lenstools_calculate_agretopicity } from '../lenstools/lenstools.nf'
include { lenstools_add_snv_metadata } from '../lenstools/lenstools.nf'

include { lenstools_make_indel_peptides } from '../lenstools/lenstools.nf'
include { lenstools_make_indel_peptides_context } from '../lenstools/lenstools.nf'
include { lenstools_add_indel_metadata } from '../lenstools/lenstools.nf'

include { fusions_to_fusion_cds_fastas } from '../seq_variation/seq_variation.nf'
include { lenstools_make_fusion_peptides } from '../lenstools/lenstools.nf'
include { lenstools_make_fusion_peptides_context } from '../lenstools/lenstools.nf'
include { lenstools_add_fusion_metadata } from '../lenstools/lenstools.nf'

include { lenstools_get_expressed_ervs_bed } from '../lenstools/lenstools.nf'
include { lenstools_get_expressed_viral_bed } from '../lenstools/lenstools.nf'
include { lenstools_get_expressed_selfs_bed } from '../lenstools/lenstools.nf'
include { bcftools_simple_germline_consensus as bcftools_simple_germline_consensus_erv } from '../bcftools/bcftools.nf'
include { bcftools_simple_germline_consensus as bcftools_simple_germline_consensus_viral } from '../bcftools/bcftools.nf'
include { bcftools_consensus_germline } from '../bcftools/bcftools.nf'
include { lenstools_make_erv_peptides } from '../lenstools/lenstools.nf'
include { lenstools_add_erv_metadata } from '../lenstools/lenstools.nf'

include { htslib_bgzip } from '../htslib/htslib.nf'
include { bcftools_index } from '../bcftools/bcftools.nf'

include { lenstools_make_self_antigen_peptides } from '../lenstools/lenstools.nf'
include { lenstools_add_self_antigen_metadata } from '../lenstools/lenstools.nf'

include { lenstools_make_viral_peptides } from '../lenstools/lenstools.nf'
include { lenstools_get_viral_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_erv_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_self_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_fusion_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_snv_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_get_indel_peptide_read_count } from '../lenstools/lenstools.nf'
include { lenstools_add_viral_metadata } from '../lenstools/lenstools.nf'

include { netmhcpan } from '../netmhcpan/netmhcpan.nf'
include { netmhcpan_rna } from '../netmhcpan/netmhcpan.nf'
include { netmhcpan as netmhcpan_wt } from '../netmhcpan/netmhcpan.nf'
include { netmhcpan as netmhcpan_mt } from '../netmhcpan/netmhcpan.nf'

include { netmhcstabpan } from '../netmhcstabpan/netmhcstabpan.nf'
include { netmhcstabpan_rna } from '../netmhcstabpan/netmhcstabpan.nf'

include { antigen_garnish_foreignness } from '../antigen.garnish/antigen.garnish.nf'
include { antigen_garnish_foreignness_rna } from '../antigen.garnish/antigen.garnish.nf'
include { antigen_garnish_dissimilarity } from '../antigen.garnish/antigen.garnish.nf'
include { antigen_garnish_dissimilarity_rna } from '../antigen.garnish/antigen.garnish.nf'

include { combine_sample_somatic_files } from '../utilities/utilities.nf'
include { combine_sample_somatic_files as combine_sample_somatic_files_consensus_alleles} from '../utilities/utilities.nf'

include { samtools_faidx_fetch } from '../samtools/samtools.nf'
include { samtools_index } from '../samtools/samtools.nf'

include { starfusion_fusions_to_extracted_reads } from '../fusion/fusion.nf'
include { seqtk_subseq } from '../seqtk/seqtk.nf'

workflow filtered_snvs_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) SNVs through NetMHCpan and other
// downstream tools.
//
// input:
// output:
//
// require:
//   VCFS_AND_INDICES
//   CLASS1_ALLELES
  take:
    vcfs
    somatic_vcfs_w_var_tx_seqs
    expressed_txs
    alleles
    pcvi_results
    quants
    rna_bams_bais
    gtf
    netmhcpan_parameters
    netmhcstabpan_parameters
    lenstools_add_snv_metadata_parameters
    ag_data_dir
    ag_species
    pep_ref
  main:
    somatic_vcfs_w_var_tx_seqs
      .join(expressed_txs, by: [0, 1, 2, 3])
      .set{ somatic_vcfs_w_var_tx_seqs_w_expressed_txs }
    lenstools_make_snv_peptides_context(
      somatic_vcfs_w_var_tx_seqs_w_expressed_txs,
      gtf,
      pep_ref)
    lenstools_make_snv_peptides_context.out.mutant_fastas
      .join(rna_bams_bais, by: [0])
      .map{ [it[0], it[1], it[2], it[5], it[3], it[4], it[8], it[9]] }
      .set{ rna_cov_inputs }
    lenstools_make_snv_peptides_context.out.wildtype_fastas
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[7]] }
      .set{ wildtype_peptides_and_alleles }
    lenstools_make_snv_peptides_context.out.mutant_fastas
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[7]] }
      .set{ mutant_peptides_and_alleles }
    netmhcpan_wt(
      wildtype_peptides_and_alleles,
      netmhcpan_parameters)
    netmhcpan_mt(
      mutant_peptides_and_alleles,
      netmhcpan_parameters)
    netmhcstabpan(
      mutant_peptides_and_alleles,
      netmhcstabpan_parameters)
    netmhcpan_wt.out.c1_antis
      .join(netmhcpan_mt.out.c1_antis, by: [0, 1, 2, 3])
      .join(lenstools_make_snv_peptides_context.out.mutant_fastas, by: [0, 1, 2, 3])
      .set{ wt_and_mt_peptides }
    lenstools_filter_mutant_snv_peptides(
        wt_and_mt_peptides,
        '11')
    lenstools_calculate_agretopicity(
        lenstools_filter_mutant_snv_peptides.out.filt_netmhcpans)
    lenstools_make_snv_peptides_context.out.nt_fastas
      .join(lenstools_calculate_agretopicity.out.agreto_netmhcpans, by: [0, 1, 2, 3])
      .set{ nt_fas_and_agretos }
    rna_bams_bais
      .join(nt_fas_and_agretos, by: 0)
      .map{ [it[0], it[5], it[6], it[1], it[2], it[3], it[4], it[8], it[9]] }
      .set{ peptide_count_inputs }
    lenstools_get_snv_peptide_read_count(
      peptide_count_inputs,
      gtf,
      params.samps_out_dir)
    antigen_garnish_foreignness(
      lenstools_get_snv_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity(
      lenstools_get_snv_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    lenstools_get_snv_peptide_read_count.out.netmhcpan_peptide_counts
      .join(pcvi_results, by: 0)
      .join(lenstools_make_snv_peptides_context.out.mutant_fastas, by: 0)
      .join(quants, by: 0)
      .join(netmhcstabpan.out.c1_stabs, by: 0)
      .join(antigen_garnish_foreignness.out.foreignness_files, by: 0)
      .join(antigen_garnish_dissimilarity.out.dissimilarity_files, by: 0)
      .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[9], it[13], it[16], it[20], it[25], it[30]] }
      .set{ metadata_snv_input }
    lenstools_add_snv_metadata(
      metadata_snv_input,
      gtf,
      lenstools_add_snv_metadata_parameters,
      params.lens_out_dir)
 emit:
   agreto_netmhcpans = lenstools_calculate_agretopicity.out.agreto_netmhcpans
   mutant_fastas = lenstools_make_snv_peptides_context.out.mutant_fastas
   metadatas = lenstools_add_snv_metadata.out.metadatas
}


workflow filtered_indels_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) InDels through netMHCpan and
// other downstream tools.
//
// input:
// output:
//
// require:
//   VCFS_AND_INDICES
//   CLASS1_ALLELES
  take:
    vcfs
    somatic_vcfs_w_var_tx_seqs
    expressed_txs
    alleles
    pcvi_results
    quants
    rna_bams_bais
    gtf
    netmhcpan_parameters
    lenstools_add_indel_metadata_parameters
    netmhcstabpan_parameters
    ag_data_dir
    ag_species
    pep_ref
  main:
    somatic_vcfs_w_var_tx_seqs
      .join(expressed_txs, by: [0, 1, 2, 3])
      .set{ somatic_vcfs_w_var_tx_seqs_w_expressed_txs }
    lenstools_make_indel_peptides_context(
      somatic_vcfs_w_var_tx_seqs_w_expressed_txs,
      gtf,
      pep_ref)
    lenstools_make_indel_peptides_context.out.peptide_fastas
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[7]] }
      .set{ peptides_and_alleles }
    netmhcpan(
      peptides_and_alleles,
      netmhcpan_parameters)
    netmhcpan.out.c1_antis
      .join(lenstools_make_indel_peptides_context.out.peptide_fastas, by: [0, 1, 2, 3])
      .set{ mt_nmp_fas }
    lenstools_filter_mutant_indel_peptides(
      mt_nmp_fas,
      '11')
    netmhcstabpan(
      peptides_and_alleles,
      netmhcstabpan_parameters)
    lenstools_make_indel_peptides_context.out.nt_fastas
//      .join(netmhcpan.out.c1_antis, by: [0, 1, 2, 3])
      .join(lenstools_filter_mutant_indel_peptides.out.filt_netmhcpans, by: [0, 1, 2, 3])
      .set{ nt_fas_and_c1_antis }
    rna_bams_bais
      .join(nt_fas_and_c1_antis, by: 0)
      .map{ [it[0], it[5], it[6], it[1], it[2], it[3], it[4], it[8], it[9]] }
      .set{ peptide_count_inputs }
    lenstools_get_indel_peptide_read_count(
      peptide_count_inputs,
      gtf,
      params.samps_out_dir)
    antigen_garnish_foreignness(
      lenstools_get_indel_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity(
      lenstools_get_indel_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    lenstools_get_indel_peptide_read_count.out.netmhcpan_peptide_counts
      .join(pcvi_results, by: 0)
      .join(lenstools_make_indel_peptides_context.out.peptide_fastas, by: 0)
      .join(quants, by: 0)
      .join(netmhcstabpan.out.c1_stabs, by: 0)
      .join(antigen_garnish_foreignness.out.foreignness_files, by: 0)
      .join(antigen_garnish_dissimilarity.out.dissimilarity_files, by: 0)
      .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[9], it[13], it[16], it[20], it[25], it[30]] }
      .set{ metadata_indel_input }
    lenstools_add_indel_metadata(
      metadata_indel_input,
      gtf,
      lenstools_add_indel_metadata_parameters,
      params.lens_out_dir)
 emit:
   c1_antis = netmhcpan.out.c1_antis
   peptide_fastas = lenstools_make_indel_peptides_context.out.peptide_fastas
   metadatas = lenstools_add_indel_metadata.out.metadatas
}


workflow filtered_fusions_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) SNVs through netMHCpan and other
// downstream tools.
//
// input:
// output:
//
// require:
//   VCFS_AND_INDICES
//   CLASS1_ALLELES
  take:
    tumor_rna_reads
    fusions
    fusion_supporting_reads
    alleles
    phased_germline_vcfs
    gtf
    dna_ref
    seqtk_subseq_parameters
    netmhcpan_parameters
    netmhcstabpan_parameters
    lenstools_get_fusion_peptide_read_count_parameters
    lenstools_add_fusion_metadata_parameters
    bcftools_index_phased_germline_parameters
    lenstools_get_fusion_transcripts_bed_parameters
    samtools_faidx_fetch_parameters
    ag_data_dir
    ag_species
  main:
    starfusion_fusions_to_extracted_reads(
      fusion_supporting_reads)
    tumor_rna_reads
      .join(starfusion_fusions_to_extracted_reads.out.fusion_read_names, by: [0, 1, 2])
      .set{ tumor_reads_and_fusion_read_names }
    seqtk_subseq(
      tumor_reads_and_fusion_read_names,
      '.fusion_reads',
      seqtk_subseq_parameters)
    fusions_to_fusion_cds_fastas(
      fusions,
      gtf,
      dna_ref,
      phased_germline_vcfs,
      bcftools_index_phased_germline_parameters,
      lenstools_get_fusion_transcripts_bed_parameters,
      samtools_faidx_fetch_parameters)
    fusions
      .join(fusions_to_fusion_cds_fastas.out.fusion_tx_exons, by:[0, 1])
      .map{ [it[0], it[1], it[4], it[2], it[3], it[6]] }.set{ fusions_w_tx_exons }
    lenstools_make_fusion_peptides_context(
      fusions_w_tx_exons,
      gtf)
    lenstools_make_fusion_peptides_context.out.fusion_peptides
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[4], it[7]] }
      .set{ mutant_peptides_and_alleles }
    netmhcpan(
      mutant_peptides_and_alleles,
      netmhcpan_parameters)
    netmhcstabpan(
      mutant_peptides_and_alleles,
      netmhcstabpan_parameters)
    fusion_supporting_reads
      .join(seqtk_subseq.out.extracted_fqs, by: [0, 1, 2])
      .join(lenstools_make_fusion_peptides_context.out.fusion_nts, by: [0, 1])
      .join(netmhcpan.out.c1_antis, by: [0, 1])
      .map{ [it[0], it[1], it[5], it[6], it[3], it[4], it[7], it[10]] }
      .set{ fusions_and_fusion_reads_and_fusion_nts_and_netmhcpan }
    lenstools_get_fusion_peptide_read_count(
      fusions_and_fusion_reads_and_fusion_nts_and_netmhcpan,
      lenstools_get_fusion_peptide_read_count_parameters)
    antigen_garnish_foreignness(
      lenstools_get_fusion_peptide_read_count.out.netmhcpan_peptide_counts.map{ [it[0], it[2], 'NA', it[1], it[3], it[4]] },
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity(
      lenstools_get_fusion_peptide_read_count.out.netmhcpan_peptide_counts.map{ [it[0], it[2], 'NA', it[1], it[3], it[4]] },
      ag_data_dir,
      ag_species)
    lenstools_get_fusion_peptide_read_count.out.netmhcpan_peptide_counts
      .join(lenstools_make_fusion_peptides_context.out.fusion_peptides, by: [0, 1, 2, 3])
      .join(fusions, by: [0, 1])
      .join(netmhcstabpan.out.c1_stabs, by: 0)
      .join(antigen_garnish_foreignness.out.foreignness_files, by: 0)
      .join(antigen_garnish_dissimilarity.out.dissimilarity_files, by: 0)
      .map{ [it[0], it[1], it[2], it[3], it[4], it[5], it[7], it[11], it[16], it[21]] }
      .set{ c1_antis_and_fusions_and_fasta }
    lenstools_add_fusion_metadata(
      c1_antis_and_fusions_and_fasta,
      lenstools_add_fusion_metadata_parameters,
      params.lens_out_dir)
 emit:
   c1_antis = netmhcpan.out.c1_antis
   fusion_peptides = lenstools_make_fusion_peptides_context.out.fusion_peptides
   metadatas = lenstools_add_fusion_metadata.out.metadatas
}


workflow filtered_ervs_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) InDels through netMHCpan and
// other downstream tools.
//
// input:
// output:
//
// require:
//   EXPRESSED_HERVS
//   CLASS1_ALLELES
  take:
    ref
    expressed_ervs
    alleles
    rna_tumor_bams
    rna_tumor_tx_bams
    rna_quants
    geve_general
    samtools_index_parameters
    lenstools_get_expressed_ervs_bed_parameters
    lenstools_make_erv_peptides_parameters
    netmhcpan_parameters
    netmhcstabpan_parameters
    lenstools_add_erv_metadata_parameters
    ag_data_dir
    ag_species
  main:
    samtools_index(
      rna_tumor_tx_bams,
      samtools_index_parameters)
    lenstools_get_expressed_ervs_bed(
      expressed_ervs,
      geve_general,
      lenstools_get_expressed_ervs_bed_parameters)
    bcftools_simple_germline_consensus_erv(
        rna_tumor_bams,
        ref,
        lenstools_get_expressed_ervs_bed.out.expressed_ervs_beds)
    expressed_ervs
      .join(bcftools_simple_germline_consensus_erv.out.consensus_fastas, by:[0, 1, 2])
      .set{ ervs_and_fastas }
    lenstools_make_erv_peptides(
      ervs_and_fastas,
      geve_general,
      lenstools_make_erv_peptides_parameters)
    lenstools_make_erv_peptides.out.erv_peptides
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[7]] }
      .set{ peptides_and_alleles }
    netmhcpan_rna(
      peptides_and_alleles,
      netmhcpan_parameters)
    netmhcstabpan_rna(
      peptides_and_alleles,
      netmhcstabpan_parameters)
    samtools_index.out.bams_and_bais
      .join(lenstools_make_erv_peptides.out.erv_nts, by: [0, 1, 2])
      .join(netmhcpan_rna.output.c1_antis, by: [0, 1, 2])
      .set{ erv_bams_and_bais_and_consensus_and_netmhcpan }
    lenstools_get_erv_peptide_read_count(
      erv_bams_and_bais_and_consensus_and_netmhcpan,
      params.samps_out_dir)
    antigen_garnish_foreignness_rna(
      lenstools_get_erv_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity_rna(
      lenstools_get_erv_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    lenstools_get_erv_peptide_read_count.out.netmhcpan_peptide_counts
      .join(lenstools_make_erv_peptides.out.erv_peptides, by: [0, 1, 2])
      .join(rna_quants, by: [0, 1, 2])
      .join(bcftools_simple_germline_consensus_erv.out.vcfs, by: [0, 1, 2])
      .join(netmhcstabpan_rna.out.c1_stabs, by: [0, 1, 2])
      .join(antigen_garnish_foreignness_rna.out.foreignness_files, by: [0, 1, 2])
      .join(antigen_garnish_dissimilarity_rna.out.dissimilarity_files, by: [0, 1, 2])
      .set{ antis_and_ervs_and_quants_and_vcfs }
    lenstools_add_erv_metadata(
      antis_and_ervs_and_quants_and_vcfs,
      geve_general,
      lenstools_add_erv_metadata_parameters,
      params.lens_out_dir)
 emit:
  c1_antis = netmhcpan_rna.out.c1_antis
   metadatas = lenstools_add_erv_metadata.out.metadatas
}


workflow filtered_selfs_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) InDels through netMHCpan and
// other downstream tools.
//
// input:
// output:
//
// require:
//   EXPRESSED_SELFS
//   CLASS1_ALLELES
  take:
    rna_tumor_tx_bams
    expressed_selfs
    alleles
    quants
    germline_vcfs
    somatic_vcfs
    gtf
    dna_ref
    cta_self_gene_list
    samtools_index_parameters
    lenstools_get_expressed_selfs_bed_parameters
    samtools_faidx_fetch_parameters
    bcftools_index_parameters
    netmhcpan_parameters
    netmhcstabpan_parameters
    ag_data_dir
    ag_species
  main:
    samtools_index(
      rna_tumor_tx_bams,
      samtools_index_parameters) // Add parameter string here.
    //This needs some attention. The run name is null currently.
    expressed_selfs.map{ [it[0], it[2], it[3]] }
      .join(germline_vcfs.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(somatic_vcfs.map{ [it[0], it[3], it[4]] }, by: [0, 1])
      .set{ selfs_with_variants }
    lenstools_get_expressed_selfs_bed(
      expressed_selfs,
      gtf,
      lenstools_get_expressed_selfs_bed_parameters)
    samtools_faidx_fetch(
      lenstools_get_expressed_selfs_bed.out.expressed_selfs_beds,
      dna_ref,
      'expressed_selfs',
      samtools_faidx_fetch_parameters)
    htslib_bgzip(
      germline_vcfs)
    bcftools_index(
      htslib_bgzip.out.bgzip_files,
      bcftools_index_parameters)
    bcftools_index.out.vcfs_w_csis
      .join(samtools_faidx_fetch.out.fetched_fastas, by: [0, 2])
      .map{ [it[0], it[5], it[1], it[3], it[4], it[6]] }
      .set{ vcfs_w_csis_w_refs }
    bcftools_consensus_germline(
      vcfs_w_csis_w_refs,
      '-H R')
    bcftools_consensus_germline.out.consensus_fastas
      .join(expressed_selfs, by: [0, 1, 2])
      .set{ consensus_fastas_and_expressed_selfs }
    lenstools_make_self_antigen_peptides(
      consensus_fastas_and_expressed_selfs,
      gtf)
    lenstools_make_self_antigen_peptides.out.self_antigen_peptides
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[6]] }
      .set{peptides_and_alleles }
    netmhcpan_rna(
      peptides_and_alleles,
      netmhcpan_parameters)
    netmhcstabpan_rna(
      peptides_and_alleles,
      netmhcpan_parameters)
    samtools_index.out.bams_and_bais
      .join(lenstools_make_self_antigen_peptides.out.self_antigen_nts, by: [0, 1, 2])
      .join(netmhcpan_rna.output.c1_antis, by: [0, 1, 2])
      .set{ self_bams_and_bais_and_consensus_and_netmhcpan }
    lenstools_get_self_peptide_read_count(
      self_bams_and_bais_and_consensus_and_netmhcpan,
      params.samps_out_dir)
    antigen_garnish_foreignness_rna(
      lenstools_get_self_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity_rna(
      lenstools_get_self_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    lenstools_get_self_peptide_read_count.out.netmhcpan_peptide_counts
      .map{ [it[0], it[2], it[3]] }
      .join(quants.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(lenstools_make_self_antigen_peptides.out.self_antigen_peptides.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(netmhcstabpan_rna.out.c1_stabs.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(antigen_garnish_foreignness_rna.out.foreignness_files.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .join(antigen_garnish_dissimilarity_rna.out.dissimilarity_files.map{ [it[0], it[2], it[3]] }, by: [0, 1])
      .map{ [it[0], it[0], it[1], it[2], it[3], it[4], it[5], it[6], it[7]] }
      .set{ c1_antis_and_quants }
    lenstools_add_self_antigen_metadata(
      c1_antis_and_quants,
      gtf,
      cta_self_gene_list,
      params.lens_out_dir)
 emit:
   c1_antis = netmhcpan_rna.out.c1_antis
   self_antigen_peptides = lenstools_make_self_antigen_peptides.out.self_antigen_peptides
   metadatas = lenstools_add_self_antigen_metadata.out.metadatas
}


workflow filtered_viruses_to_class1_neos {
// Runs filtered (e.g. expressed and isolated) InDels through netMHCpan and
// other downstream tools.
//
// input:
// output:
//
// require:
//   VCFS_AND_INDICES
//   CLASS1_ALLELES
  take:
    viruses
    virus_bams_and_bais
    alleles
    gff
    cds_ref
    pep_ref
    lenstools_get_expressed_viral_bed_parameters
    lenstools_make_viral_peptides_parameters
    netmhcpan_parameters
    netmhcstabpan_parameters
    lenstools_add_viral_metadata_parameters
    ag_data_dir
    ag_species
  main:
    lenstools_get_expressed_viral_bed(
      viruses,
      cds_ref,
      lenstools_get_expressed_viral_bed_parameters)
    bcftools_simple_germline_consensus_viral(
      virus_bams_and_bais,
      cds_ref,
      lenstools_get_expressed_viral_bed.out.expressed_viral_beds)
    viruses
      .join(virus_bams_and_bais, by: [0, 1, 2])
      .set{ viruses_and_bams_and_bais }
    lenstools_make_viral_peptides(
      bcftools_simple_germline_consensus_viral.out.consensus_fastas,
      lenstools_make_viral_peptides_parameters)
    lenstools_make_viral_peptides.out.viral_peptides
      .combine(alleles, by:[0])
      .map{ [it[0], it[1], it[2], it[3], it[6]] }
      .set{ peptides_and_alleles }
    netmhcpan_rna(
      peptides_and_alleles,
      netmhcpan_parameters)
    netmhcstabpan_rna(
      peptides_and_alleles,
      netmhcpan_parameters)
    virus_bams_and_bais
      .join(bcftools_simple_germline_consensus_viral.out.consensus_fastas, by: [0, 1, 2])
      .join(netmhcpan_rna.output.c1_antis, by: [0, 1, 2])
      .set{ virus_bams_and_bais_and_consensus_and_netmhcpan }
    lenstools_get_viral_peptide_read_count(
      virus_bams_and_bais_and_consensus_and_netmhcpan,
      params.samps_out_dir)
    antigen_garnish_foreignness_rna(
      lenstools_get_viral_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    antigen_garnish_dissimilarity_rna(
      lenstools_get_viral_peptide_read_count.out.netmhcpan_peptide_counts,
      ag_data_dir,
      ag_species)
    lenstools_get_viral_peptide_read_count.out.netmhcpan_peptide_counts
      .join(netmhcstabpan_rna.out.c1_stabs, by:[0, 1, 2])
      .join(antigen_garnish_foreignness_rna.out.foreignness_files, by: [0, 1, 2])
      .join(antigen_garnish_dissimilarity_rna.out.dissimilarity_files, by: [0, 1, 2])
      .set{ all_viral_metadata_sources }
    lenstools_add_viral_metadata(
      all_viral_metadata_sources,
      cds_ref,
      pep_ref,
      lenstools_add_viral_metadata_parameters,
      params.lens_out_dir)
 emit:
   c1_antis = netmhcpan_rna.out.c1_antis
   viral_peptides = lenstools_make_viral_peptides.out.viral_peptides
   metadatas = lenstools_add_viral_metadata.out.metadatas
}
