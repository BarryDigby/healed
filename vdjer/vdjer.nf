#!/usr/bin/env nextflow

process vdjer {
// Runs vdjer
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(prefix) - FASTQ Prefix
//     val(dataset) - Dataset
//     path(aln) - Alignment File
//     path(bai) - Alignment Index File
//   path(vdjer_refs) - V'DJer Reference Directory
//   val(ins) - Insert Size
//   val(chain) - Chain (IGH, IGL, IGV)
//   val(out_dir) - Output Directory
//
// output:
//

  label 'vdjer'
  publishDir "${out_dir}/${dataset}/${pat_name}/${prefix}", pattern: "vdjer-*", mode: 'copyNoFollow'
  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(aln), path(bai)
  path vdjer_refs
  val ins
  val chain
  val out_dir

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*'), emit: vdj_calls
  path('vdjer-*')

  script:
  """
  LOWER=\$(echo ${chain} | tr '[:upper:]' '[:lower:]')
  UPPER=\$(echo ${chain} | tr '[:lower:]' '[:upper:]')

  CMD="vdjer --in ${aln} --t ${task.cpus} --ins ${ins} --chain \${UPPER} --ref-dir ${vdjer_refs}/\${LOWER}"
  eval \${CMD}

  #For publishing output
  if [[ ! -z ${out_dir} ]]; then
    VERSION=\$(echo ${task.container} | rev | cut -f 1 -d ':' | rev)
    CMD_MD5=\$(echo \${CMD} | md5sum | cut -f 1 -d ' ' | cut -c1-5)
    ALN_MD5=\$(md5sum ${aln} | cut -f 1 -d ' ' | cut -c1-5)
    OUTDIR="vdjer-\${VERSION}=\${CMD_MD5}/${aln}=\${ALN_MD5}"
    mkdir -p \${OUTDIR}
    ln -s \${PWD}/* \${OUTDIR}/
    cp .command.* \${OUTDIR}/
  fi
  """
}

workflow vdjer_all {
  take:
    refs
    alns
    ins
  main:
    vdjer_igh(refs, alns, ins)
    vdjer_igl(refs, alns, ins)
    vdjer_igk(refs, alns, ins)
}

workflow vdjer_igh {
  take:
    refs
    alns
    ins
  main:
    vdjer(refs, alns, ins, 'igh', params.samps_out_dir)
  emit:
    vdjer.out.vdj_calls
}

workflow vdjer_igl {
  take:
    refs
    alns
    ins
  main:
    vdjer(refs, alns, ins, 'igl', params.samps_out_dir)
  emit:
    vdjer.out.vdj_calls
}

workflow vdjer_igk {
  take:
    refs
    alns
    ins
  main:
    vdjer(refs, alns, ins, 'igk', params.samps_out_dir)
  emit:
    vdjer.out.vdj_calls
}
