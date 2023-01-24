#!/usr/bin/env nextflow

include { hervquant_filter_alignments } from '../hervquant/hervquant.nf'
include { manifest_to_rna_procd_fqs } from '../preproc/preproc.nf'
include { salmon_aln_quant } from '../salmon/salmon.nf'
include { star_index } from '../star/star.nf'
include { star_index as star_index_human } from '../star/star.nf'
include { star_index as star_index_virus } from '../star/star.nf'
include { star_map } from '../star/star.nf'
include { star_map as star_map_human } from '../star/star.nf'
include { star_map as star_map_virus } from '../star/star.nf'
include { awk_get_unaligned } from '../virdetect/virdetect.nf'
include { awk_unaligned_bam_to_fastqs } from '../virdetect/virdetect.nf'
include { count_star_viral_alignments } from '../virdetect/virdetect.nf'

workflow manifest_to_hervquant {
// require:
//   MANIFEST
  take:
    manifest
  main:
    manifest_to_rna_procd_fqs(manifest)
    procd_fqs_to_hervquant(manifest_to_rna_procd_fqs.out.procd_fqs)
  emit:
    herv_counts = procd_fqs_to_hervquant.out.herv_counts
}


workflow procd_fqs_to_hervquant {
// require:
//   FQS
  take:
    fqs
  main:
    star_index(
      params.viral$procd_fqs_to_hervquant$hervquant_star_ref,
      "--limitGenomeGenerateRAM 53667544448 --genomeSAindexNbases 7")
    star_map(
      fqs,
      star_index.out.idx_files,
      "--outFilterMultimapNmax 10 --outFilterMismatchNmax 7",
      '',
      params.dummy_file)
    hervquant_filter_alignments(
      star_map.out.alns)
    salmon_aln_quant(
      hervquant_filter_alignments.out.bams,
      params.viral$procd_fqs_to_hervquant$hervquant_salmon_ref,
      '')
  emit:
    herv_counts = salmon_aln_quant.out.quants
    herv_bams = hervquant_filter_alignments.out.bams
}


workflow manifest_to_virdetect {
// require:
//   MANIFEST
  take:
    manifest
  main:
    manifest_to_rna_procd_fqs(manifest)
    procd_fqs_to_virdetect(manifest_to_rna_procd_fqs.out.procd_fqs)
//  emit:
}


workflow procd_fqs_to_virdetect {
// require:
//   params.viral$procd_fqs_to_virdetect$fqs
  take:
    fqs
  main:
    star_index_human(
      params.viral$procd_fqs_to_virdetect$virdetect_host_ref,
      "--genomeChrBinNbits 14 --sjdbGTFfile ${params.viral$procd_fqs_to_virdetect$gtf}")
    star_index_virus(
      params.viral$procd_fqs_to_virdetect$virdetect_viral_ref,
      '--genomeSAindexNbases 7')
    star_map_human(
      fqs,
      star_index_human.out.idx_files,
      '--outFilterMultimapNmax 1000 --outSAMunmapped Within --limitOutSAMoneReadBytes 1000000',
      '',
      params.dummy_file)
    awk_get_unaligned(star_map_human.out.alns)
    awk_unaligned_bam_to_fastqs(awk_get_unaligned.out.bams)
    star_map_virus(
      awk_unaligned_bam_to_fastqs.out.fqs,
      star_index_virus.out.idx_files,
      '--outFilterMismatchNmax 4 --outFilterMultimapNmax 1000 --limitOutSAMoneReadBytes 1000000 --outSAMtype BAM SortedByCoordinate',
      '',
      params.dummy_file)
    count_star_viral_alignments(star_map_virus.out.alns)
  emit:
    viral_counts = count_star_viral_alignments.out.viral_counts
    unaligned_fqs =  awk_unaligned_bam_to_fastqs.out.fqs
}


workflow alns_to_virdetect {
// require:
//   params.viral$alns_to_virdetect$alns
  take:
    alns
    virdetect_viral_ref
  main:
    star_index_virus(
      virdetect_viral_ref,
      '--genomeSAindexNbases 7')
    awk_get_unaligned(alns)
    awk_unaligned_bam_to_fastqs(awk_get_unaligned.out.bams)
    star_map_virus(
      awk_unaligned_bam_to_fastqs.out.fqs,
      star_index_virus.out.idx_files,
      '--outFilterMismatchNmax 4 --outFilterMultimapNmax 1000 --limitOutSAMoneReadBytes 1000000 --outSAMtype BAM SortedByCoordinate',
      '',
      params.dummy_file)
    count_star_viral_alignments(star_map_virus.out.alns)
  emit:
    viral_alns = star_map_virus.out.alns
    viral_counts = count_star_viral_alignments.out.viral_counts
    unaligned_fqs =  awk_unaligned_bam_to_fastqs.out.fqs
}

workflow unaligned_fqs_to_virdetect_cds_counts {
// require:
//   UNALIGNED_READS
  take:
    fqs
    viral_cds_ref
  main:
    star_index_virus(
      viral_cds_ref,
      '--genomeSAindexNbases 7')
    star_map_virus(
      fqs,
      star_index_virus.out.idx_files,
//      '--outFilterMismatchNmax 4 --outFilterMultimapNmax 1000 --limitOutSAMoneReadBytes 1000000 --outSAMtype BAM SortedByCoordinate',
      '--outFilterMismatchNmax 3 --outFilterMultimapNmax 1 --limitOutSAMoneReadBytes 1000000 --outSAMtype BAM SortedByCoordinate',
      '',
      params.dummy_file)
    count_star_viral_alignments(star_map_virus.out.alns)
  emit:
    viral_cds_alns = star_map_virus.out.alns
    viral_cds_counts = count_star_viral_alignments.out.viral_counts
}
