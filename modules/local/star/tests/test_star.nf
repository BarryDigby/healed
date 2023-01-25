nextflow.enable.dsl=2

include { star_index } from '../star.nf'
include { star_map } from '../star.nf'

params.pwd = ''

workflow test_star_index {
  main:
    star_index(
      "${params.pwd}/data/test_star_index/NC_045512.2.fa",
      '')
}

workflow test_star_map {
  main:
    //star_index is being included as part of this workflow since it outputs
    //are large (~1.5G) even with the relatively small reference file.
    star_index(
      "${params.pwd}/data/test_star_index/NC_045512.2.fa",
      '')
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_star_map/SRR18009685_1.fastq.gz",
                "${params.pwd}/data/test_star_map/SRR18009685_2.fastq.gz"])
           .set{ sample_inputs }
    star_map(
      sample_inputs,
      star_index.out.idx_files,
      '--quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --twopassMode Basic',
      'toTranscriptome',
      "${params.pwd}/data/test_star_map/GCF_009858895.2_ASM985889v3_genomic.gtf")
}
