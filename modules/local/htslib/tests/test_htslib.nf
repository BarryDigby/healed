nextflow.enable.dsl=2

@Grab(group='commons-io', module='commons-io', version='2.11.0')

import org.apache.commons.codec.digest.DigestUtils

include { htslib_bgzip_somatic } from '../htslib.nf'
include { htslib_bgzip } from '../htslib.nf'
include { htslib_tabix } from '../htslib.nf'
include { htslib_bgzip_ref } from '../htslib.nf'

params.pwd = ''

workflow test_htslib_bgzip_somatic {
  main:
    Channel.of(['test_patient',
                'test_norm_run',
                'test_tum_run',
                'test_dataset',
                "${params.pwd}/data/test_htslib_bgzip_somatic/test.vcf"])
           .set{ test_htslib_bgzip_somatic_inputs }
    htslib_bgzip_somatic(
      test_htslib_bgzip_somatic_inputs)
}

workflow test_htslib_bgzip {
  main:
    Channel.of(['test_patient',
                'test_run',
                'test_dataset',
                "${params.pwd}/data/test_htslib_bgzip/test.vcf"])
           .set{ test_htslib_bgzip_inputs }
    htslib_bgzip(
      test_htslib_bgzip_inputs)
}

workflow test_htslib_tabix {
  main:
    Channel.of("${params.pwd}/data/test_htslib_tabix/test.vcf.gz")
           .set{ test_htslib_tabix_inputs }
    htslib_tabix(
      test_htslib_tabix_inputs)
}

workflow {
  test_htslib_bgzip_somatic()
  test_htslib_bgzip()
  test_htslib_tabix()
}
