nextflow.enable.dsl=2

include { deepvariant } from '../deepvariant.nf'

params.pwd = ''

workflow test_deepvariant {
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_deepvariant/test.bam",
                "${params.pwd}/data/test_deepvariant/test.bam.bai"])
           .set{ sample_inputs }
    Channel.of(["${params.pwd}/data/test_deepvariant/chr1.fa",
                "${params.pwd}/data/test_deepvariant/chr1.fa.fai",
                "${params.pwd}/data/test_deepvariant/dummy_file"])
           .set{ ref_inputs }
    deepvariant(
      sample_inputs,
      ref_inputs,
      "${params.pwd}/data/test_deepvariant/test.bed",
      "WES",
      '',
      ".deepv")
}
 
      



