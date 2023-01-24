nextflow.enable.dsl=2

include { bcftools_sort } from '../bcftools.nf'
include { bcftools_consensus } from '../bcftools.nf'
include { bcftools_filter } from '../bcftools.nf'
include { bcftools_index } from '../bcftools.nf'
include { bcftools_index_somatic } from '../bcftools.nf'
include { bcftools_stats } from '../bcftools.nf'
include { bcftools_stats_somatic } from '../bcftools.nf'

params.pwd = ''

workflow test_bcftools_sort {
  main:
    Channel.of(['test_patient',
                'test_run', 
                'test_dataset', 
                "${params.pwd}/data/test_bcftools_sort/test_1.vcf"])
           .set{ bcftools_sort_inputs }
    bcftools_sort(
      bcftools_sort_inputs,
      '')
}

workflow test_bcftools_sort_gzipped {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_sort/test_2.vcf.gz"])
           .set{ bcftools_sort_inputs }
    bcftools_sort(
      bcftools_sort_inputs,
      '')
}

workflow test_bcftools_consensus_snv {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tum_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_consensus/test_bcftools_consensus_snv/test.vcf.gz",
                "${params.pwd}/data/test_bcftools_consensus/test_bcftools_consensus_snv/test.vcf.gz.csi",
                "${params.pwd}/data/test_bcftools_consensus/test_bcftools_consensus_snv/test.fa"])
           .set{ bcftools_consensus_inputs }
    bcftools_consensus(
      bcftools_consensus_inputs,
      'test',
      '')
}

workflow test_bcftools_filter {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tum_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_filter/test.vcf"])
           .set{ bcftools_filter_inputs }
    bcftools_filter(
      bcftools_filter_inputs,
      '-i \'FILTER=\"PASS\"\'')
}

workflow test_bcftools_index {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_index/test.vcf.gz"])
           .set{ bcftools_index_inputs }
    bcftools_index(
      bcftools_index_inputs,
      '')
}

workflow test_bcftools_index_somatic {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tum_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_index_somatic/test.vcf.gz"])
           .set{ bcftools_index_somatic_inputs }
    bcftools_index_somatic(
      bcftools_index_somatic_inputs,
      '')
}

workflow test_bcftools_stats {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_stats/test.vcf.gz"])
           .set{ bcftools_stats_inputs }
    bcftools_stats(
      bcftools_stats_inputs,
      '')
}

workflow test_bcftools_stats_somatic {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tum_run',
                'test_dataset',
                "${params.pwd}/data/test_bcftools_stats_somatic/test.vcf.gz"])
           .set{ bcftools_stats_somatic_inputs }
    bcftools_stats_somatic(
      bcftools_stats_somatic_inputs,
      '')
}


workflow {
  test_bcftools_sort()
  test_bcftools_sort_gzipped()
  test_bcftools_consensus_snv()
  test_bcftools_filter()
  test_bcftools_index()
  test_bcftools_index_somatic()
  test_bcftools_stats()
  test_bcftools_stats_somatic()
}
