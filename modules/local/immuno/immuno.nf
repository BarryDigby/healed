#!/usr/bin/env nextflow

include { hlaprofiler_predict } from '../hlaprofiler/hlaprofiler.nf'

include { optitype } from '../optitype/optitype.nf'

include { mixcr_align } from '../mixcr/mixcr.nf'
include { mixcr_assemble } from '../mixcr/mixcr.nf'
include { mixcr_export } from '../mixcr/mixcr.nf'
include { mixcr_shotgun } from '../mixcr/mixcr.nf'

include { procd_fqs_to_star_alns } from '../rna_quant/rna_quant.nf'

include { samtools_index } from '../samtools/samtools.nf'

include { seqtk_sample } from '../seqtk/seqtk.nf'
include { trim_galore } from '../trim_galore/trim_galore.nf'
include { trim_galore_hlap } from '../trim_galore/trim_galore.nf'

include { get_fastqs } from '../utilities/utilities.nf'

include { vdjer } from '../vdjer/vdjer.nf'
include { vdjer_all } from '../vdjer/vdjer.nf'


workflow manifest_to_hlaprofiler_calls {
// require:
//   MANIFEST
//   params.immuno$raw_fqs_to_hlaprofiler_calls$trim_galore_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_read_count
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_seed
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_sample_suffix
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$hlaprofiler_parameters
  take:
    manifest
    trim_galore_parameters
    seqtk_read_count
    seqtk_seed
    seqtk_sample_suffix
    seqtk_parameters
    trim_galore_hlap_parameters
    hlaprofiler_parameters
  main:
    get_fastqs(
      manifest,
      params$fq_dir)
    raw_fqs_to_hlaprofiler_calls(
      get_fastqs.out.fastqs)
  emit:
    calls = raw_fqs_to_hlaprofiler_calls.out.calls
}


workflow raw_fqs_to_hlaprofiler_calls {
// require:
//   FQS
//   params.immuno$raw_fqs_to_hlaprofiler_calls$trim_galore_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_read_count
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_seed
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_sample_suffix
//   params.immuno$raw_fqs_to_hlaprofiler_calls$seqtk_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters
//   params.immuno$raw_fqs_to_hlaprofiler_calls$hlaprofiler_parameters
  take:
    fqs
    trim_galore_parameters
    seqtk_read_count
    seqtk_seed
    seqtk_sample_suffix
    seqtk_parameters
    trim_galore_hlap_parameters
    hlaprofiler_parameters
  main:
    trim_galore(
      fqs,
      params.immuno$raw_fqs_to_hlaprofiler_calls$trim_galore_parameters)
    procd_fqs_to_hlaprofiler_calls(
      trim_galore.out.procd_fqs,
      seqtk_read_count,
      seqtk_seed,
      seqtk_sample_suffix,
      seqtk_parameters,
      trim_galore_hlap_parameters,
      hlaprofiler_parameters)
  emit:
    calls = procd_fqs_to_hlaprofiler_calls.out.calls
}


workflow procd_fqs_to_hlaprofiler_calls {
// require:
//   FQS
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_read_count
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_seed
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_sample_suffix
//   params.immuno$procd_fqs_to_hlaprofiler_calls$seqtk_parameters
//   params.immuno$procd_fqs_to_hlaprofiler_calls$trim_galore_hlap_parameters
//   params.immuno$procd_fqs_to_hlaprofiler_calls$hlaprofiler_parameters
  take:
    fqs
    seqtk_read_count
    seqtk_seed
    seqtk_sample_suffix
    seqtk_parameters
    trim_galore_hlap_parameters
    hlaprofiler_parameters
  main:
    seqtk_sample(
      fqs,
      seqtk_read_count,
      seqtk_seed,
      seqtk_sample_suffix,
      seqtk_parameters)
    trim_galore_hlap(
      seqtk_sample.out.subd_fqs,
      trim_galore_hlap_parameters) 
    hlaprofiler_predict(
      trim_galore_hlap.out.procd_fqs,
      hlaprofiler_parameters)
  emit:
    calls = hlaprofiler_predict.out.calls
}


workflow manifest_to_mixcr_clonotypes {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest,
      params.immuno$fq_dir)
    raw_fqs_to_mixcr_clonotypes(
      get_fastqs.out.fastqs)
  emit:
    clonotypes = raw_fqs_to_mixcr_clonotypes.out.clonotypes
}


workflow raw_fqs_to_mixcr_clonotypes {
// require:
//   FQS
  take:
    fqs
  main:
    trim_galore(
      fqs,
      params.immuno$raw_fqs_to_mixcr_clonotypes$trim_galore_parameters)
    procd_fqs_to_mixcr_clonotypes(
      trim_galore.out.procd_fqs)
  emit:
    clonotypes = procd_fqs_to_mixcr_clonotypes.out.clonotypes
}


workflow procd_fqs_to_mixcr_clonotypes {
// require:
//   FQS
  take:
    fqs
  main:
    mixcr_align(
      fqs,
      params.immuno$procd_fqs_to_mixcr_clonotypes$mixcr_align_parameters)
    vdjcas_to_mixcr_clonotypes(
      mixcr_align.out.vdjca_files)
  emit:
    clonotypes = vdjcas_to_mixcr_clonotypes.out.clonotypes
}


workflow vdjcas_to_mixcr_clonotypes {
// require:
//   VDJCAS
  take:
    vdjcas
  main:
    mixcr_assemble(
      vdjcas,
      params.immuno$alns_to_mixcr_clonotypes$mixcr_assemble_if_config,
      params.immuno$alns_to_mixcr_clonotypes$mixcr_assemble_parameters)
    clns_to_mixcr_clonotypes(
      mixcr_assemble.out.cln_files)
  emit:
    clonotypes = clns_to_mixcr_clonotypes.out.clonotypes
}


workflow clns_to_mixcr_clonotypes {
// require:
//   params.immuno$clns_to_mixcr_clonotypes$clns
  take:
    cln_files
  main:
    mixcr_export(
      cln_files,
      params.immuno$clns_to_mixcr_clonotypes$exports,
      params.immuno$clns_to_mixcr_clonotypes$mixcr_export_parameters)
  emit:
    clonotypes = mixcr_export.out.exported_files
}


workflow manifest_to_vdjer_calls {
// require:
//   MANIFEST
//   params.immuno$manifest_to_vdjer_calls$genome_ref
//   params.immuno$manifest_to_vdjer_calls$vdj_ref
//   params.immuno$manifest_to_vdjer_calls$chain
  take:
    manifest
    g_ref
    vdj_ref
    chain
  main:
    get_fastqs(
      manifest,
      params.immuno$fq_dir)
    raw_fqs_to_vdjer_calls(
      get_fastqs.out.fastqs,
      g_ref,
      vdj_ref,
      chain)
  emit:
    vdj_calls = raw_fqs_to_vdjer_calls.out.vdj_calls
}


workflow raw_fqs_to_vdjer_calls {
// require:
//   FQS
//   params.immuno$raw_fqs_to_vdjer_calls$genome_ref
//   params.immuno$raw_fqs_to_vdjer_calls$vdj_ref
//   params.immuno$raw_fqs_to_vdjer_calls$chain
  take:
    fqs
    g_ref
    vdj_ref
    chain
  main:
    trim_galore(
      fqs,
      params.immuno$raw_fqs_to_vdjer_calls$trim_galore_parameters)
    procd_fqs_to_vdjer_calls(
      trim_galore.out.procd_fqs,
      g_ref,
      vdj_ref,
      chain)
  emit:
    vdj_calls = procd_fqs_to_vdjer_calls.out.vdj_calls
}


workflow procd_fqs_to_vdjer_calls {
// require:
//   FQS
//   params.immuno$procd_fqs_to_vdjer_calls$genome_ref
//   params.immuno$procd_fqs_to_vdjer_calls$vdj_ref
//   params.immuno$procd_fqs_to_vdjer_calls$chain
  take:
    fqs
    g_ref
    vdj_ref
    chain
  main:
    procd_fqs_to_star_alns(
      g_ref,
      fqs)
    alns_to_vdjer_calls(
      procd_fqs_to_star_alns.out.alns,
      vdj_ref,
      chain)
  emit:
    vdj_calls = alns_to_vdjer_calls.out.vdj_calls
}


workflow alns_to_vdjer_calls {
// require:
//   ALNS
//   params.immuno$alns_to_vdjer_calls$vdj_ref
//   params.immuno$alns_to_vdjer_calls$chain
  take:
    alns
    vdj_ref
    chain
  main:
    samtools_index(alns)
    alns.join(samtools_index.out.bais, by: [0, 1, 2]).set{ alns_w_idxs }
    Channel.empty().set{ calls }
    if(chain == 'all') {
      vdjer_all(
        alns_w_idxs,
        vdj_ref,
        params.immuno$alns_to_vdjer_calls$vdjer_insert_size)
      vdjer_all.out.vdj_calls
    } else {
      vdjer(
        alns_w_idxs,
        vdj_ref,
        params.immuno$alns_to_vdjer_calls$vdjer_insert_size,
        chain)
      vdjer.out.vdj_calls.set{ calls }
    }
  emit:
    vdj_calls = calls
}


workflow manifest_to_mixcr_shotgun_clonotypes {
// require:
//   MANIFEST
  take:
    manifest
  main:
    get_fastqs(
      manifest,
      params.immuno$fq_dir)
    raw_fqs_to_mixcr_shotgun_clonotypes(
      get_fastqs.out.fastqs)
  emit:
    clonotypes = raw_fqs_to_mixcr_shotgun_clonotypes.out.clonotypes
    report = raw_fqs_to_mixcr_shotgun_clonotypes.out.report
    aligned_r1 = raw_fqs_to_mixcr_shotgun_clonotypes.out.aligned_r1
    aligned_r2 = raw_fqs_to_mixcr_shotgun_clonotypes.out.aligned_r2
}


workflow raw_fqs_to_mixcr_shotgun_clonotypes {
// require:
//   FQS
  take:
    fqs
  main:
    trim_galore(
      fqs,
      params.immuno$raw_fqs_to_mixcr_shotgun_clonotypes$trim_galore_parameters)

    procd_fqs_to_mixcr_shotgun_clonotypes(
      trim_galore.out.procd_fqs)
  emit:
    clonotypes = procd_fqs_to_mixcr_shotgun_clonotypes.out.clonotypes
    report = procd_fqs_to_mixcr_shotgun_clonotypes.out.report
    aligned_r1 = procd_fqs_to_mixcr_shotgun_clonotypes.out.aligned_r1
    aligned_r2 = procd_fqs_to_mixcr_shotgun_clonotypes.out.aligned_r2
}


workflow procd_fqs_to_mixcr_shotgun_clonotypes {
// require:
//   FQS
  take:
    fqs
  main:
    mixcr_shotgun(
      fqs,
      params.immuno$procd_fqs_to_mixcr_shotgun_clonotypes$mixcr_shotgun_parameters,
      params.immuno$procd_fqs_to_mixcr_shotgun_clonotypes$regex)
  emit:
    clonotypes = mixcr_shotgun.out.clonotypes
    report = mixcr_shotgun.out.report
    aligned_r1 = mixcr_shotgun.out.aligned_r1
    aligned_r2 = mixcr_shotgun.out.aligned_r2
}


process hlaprofiler_alleles_to_netmhcpan_alleles {
// require:
//   immuno$hlaprofiler_alleles_to_netmhcpan_alleles$hlaprofiler_calls

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(hlaprofiler_calls)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.netmhcpan.alleles"), emit: netmhcpan_alleles

  script:
  """
  grep '\\s[ABC]\\*' ${hlaprofiler_calls} |\
  cut -f 3,4 | sed 's/\\t/\\n/g' |\
  cut -f 1,2 -d ':' | sed 's/^/HLA-/' |\
  tr -d '\\*' | grep -v 'N' | sort | uniq | \
  sed -z 's/\\n/,/g' | sed 's/,\$//' |\
  sed 's/_updated//g' | sed 's/_novel//g'\
  > ${dataset}-${pat_name}-${run}.netmhcpan.alleles
  """
}

workflow extract_alleles_from_manifest {
// Extract alleles from from a manifest channel. This is useful when MHC
// alleles are known or when working with non-human species.
//
// require:
//   ALLELES
  take:
    alleles
  main:
    Channel.fromPath(alleles).splitCsv(header: true, sep: '\t')
    .map{ row -> tuple("${row.Patient_Name}", "${row.File_Prefix}", "${row.Dataset}", "${row.Alleles}") }
    .set{ patient_alleles }
  emit:
    patient_alleles
}

process user_provided_alleles_to_netmhcpan_alleles {

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), val(alleles)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.netmhcpan.alleles"), emit: netmhcpan_alleles

  script:
  """
  echo ${alleles} > ${dataset}-${pat_name}-${run}.netmhcpan.alleles
  """
}

process optitype_alleles_to_netmhcpan_alleles {
// require:
//   immuno$optitype_alleles_to_netmhcpan_alleles$optitype_calls

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(optitype_calls)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*.netmhcpan.alleles"), emit: netmhcpan_alleles

  script:
  """
  grep ^0 ${optitype_calls} |\
  grep ^0 |\
  cut -f 2,3,4,5,6,7 |\
  sed 's/^/HLA-/g' |\
  sed 's/\t/,HLA-/g' |\
  tr -d '\\*'  > ${dataset}-${pat_name}-${run}.netmhcpan.alleles
  """
}

workflow procd_fqs_to_optitype_calls {
// require:
//   FQS
  take:
    fqs
  main:
    optitype(
      fqs)
  emit:
    calls = optitype.out.calls
}
