#!/usr/bin/env nextflow

include { gtf_to_genepred }from '../gtftogenepred/gtftogenepred.nf'
include { samtools_faidx } from '../samtools/samtools.nf'

workflow gtf_to_refflat {
// Converts a GTF file to refFlat file.

// require:
//   GTF
  take:
    gtf
  main:
    gtf_to_genepred(
      gtf,
      '-genePredExt -geneNameAsName2')
    genepred_to_refflat(
        gtf_to_genepred.out.genepred)
  emit:
    refflat = genepred_to_refflat.out.refflat
}


workflow gtf_to_rrna_intervals {
// Converts a GTF file to rRNA intervals file.
//

// require:
//   FASTA
//   GTF
  take:
    fa
    gtf
    samtools_faidx_parameters
  main:
    get_chrom_sizes(
      fa,
      samtools_faidx_parameters)
    gtf_to_rrna_intervals_sub(
      gtf,
      get_chrom_sizes.out.chrom_sizes)
  emit:
    rrna_intervals = gtf_to_rrna_intervals_sub.out.rrna_intervals
}


process genepred_to_refflat {
// require:
//   GENEPRED

  tag "${genepred}"

  input:
  path genepred

  output:
  path "refFlat.txt", emit: refflat

  script:
  """
  paste <(cut -f 12 ${genepred}) <(cut -f 1-10 ${genepred}) > refFlat.txt
  """
}


process gtf_to_rrna_intervals_sub {
// require:
//   GTF
//   CHROM_SIZES

  tag "${gtf}"

  input:
  path gtf
  path chrom_sizes

  output:
  path "rrna_intervals.txt", emit: rrna_intervals

  script:
  """
  #From https://gist.github.com/ag1805x/55ba0d88f317f63423a4e54adc46eb1f
  # Sequence names and lengths. (Must be tab-delimited.)
  # Intervals for rRNA transcripts.
  perl -lane 'print "\\@SQ\\tSN:\$F[0]\\tLN:\$F[1]\\tAS:GRCh38"' ${chrom_sizes} >> rrna_intervals.txt; grep -E '(gene_biotype "rRNA"|gene_type "rRNA")' ${gtf} | awk '\$3 == "gene"' | cut -f1,4,5,7,9 | perl -lane '/gene_id "([^"]+)"/ or die "no gene_id on \$."; print join "\\t", (@F[0,1,2,3], \\\$1)' | sort -k1V -k2n -k3n >> rrna_intervals.txt
  """
}



process get_chrom_sizes_sub {
// require:
//   FASTA

  tag "${fai}"

  input:
  tuple path(fa), path(fai)

  output:
  path "chrom.sizes", emit: chrom_sizes

  script:
  """
  cut -f 1,2 ${fai} > chrom.sizes
  """
}

workflow get_chrom_sizes {
// require:
//   FASTA
//   params.ref_utiles$get_chrom_sizes$samtools_faidx_parameters
  take:
    fa
    samtools_faidx_parameters
  main:
    samtools_faidx(
      fa,
      samtools_faidx_parameters)
    get_chrom_sizes_sub(
      samtools_faidx.out.faidx_file)
  emit:
    chrom_sizes = get_chrom_sizes_sub.out.chrom_sizes
}


process get_callable_bed {
// require:
//   FASTA

  tag "${fa}"

  input:
  path(fa)

  output:
  path "callable.bed", emit: bed

  script:
  """
  grep \\> ${fa} | grep chromosome | cut -f 4,5,6 -d ':' | sed 's/:/\t/'g > callable.bed
  """
}



process extract_chroms_from_gff {
// require:
//   GFF

  tag "${gff}"

  input:
  path(gff)

  output:
  path "chroms.list", emit: chroms_list

  script:
  """
  cut -f 1 ${gff} | grep -v ^# | sort | uniq > chroms.list
  """
}



process extract_chroms_from_bed {
// require:
//   BED

  tag "${bed}"

  input:
  path bed

  output:
  path "chroms.list", emit: chroms_list

  script:
  """
  cut -f 1 ${bed} | sort | uniq > chroms.list
  """
}
