#!/usr/bin/env nextflow

process kallisto_index {
// Runs kallisto index
//
// input:
//   path fa - Reference FASTA
//   val params - Additional Parameters
//   val out_dir - Output Directory
//   val shared_dir - Shared Output Directory
//
// output:
//   tuple => emit: idx_files
//     path(fa) - Reference FASTA
//     path("${fa}*index") - Index Files
//   path('kallisto-*') - For publishDir

  label 'kallisto'
  publishDir "${out_dir}/${fa}", pattern: "kallisto-*", mode: 'copyNoFollow'
  publishDir "${shared_dir}/${fa}", pattern: "kallisto-*", mode: 'copyNoFollow'
  tag "${fa}"

  input:
  path fa
  val params
  val out_dir
  val shared_dir

  output:
  tuple path(fa), path("${fa}*index"), emit: idx_files
  path('kallisto-*')

  script:
  """
  CMD="kallisto index -i ${fa}.index ${params} ${fa}"

  VERSION=\$(kallisto | grep ^kallisto | cut -f 2 -d ' ')
  CMD_MD5=\$(echo \${CMD} | md5sum | cut -f 1 -d ' ' | cut -c1-5)
  FA_MD5=\$(md5sum ${fa} | cut -f 1 -d ' ' | cut -c1-5)

  OUTDIR="kallisto-\${VERSION}=\${CMD_MD5}/${fa}=\${FA_MD5}"
  EXP_SHARED_DIR="${shared_dir}/${out_dir}/${fa}/\${OUTDIR}"

  mkdir -p \${OUTDIR}

  if [[ ! -d \${EXP_SHARED_DIR} ]]; then
    touch .real
    eval \${CMD}
    ln -s \${PWD}/*index \${OUTDIR}/
    ln -s .command.* \${OUTDIR}/
  else
    touch .copy
    cp -r \${EXP_SHARED_DIR}/* \${OUTDIR}
    cp -r \${EXP_SHARED_DIR}/* .
  fi
  """
}


process kallisto_quant {
// Runs kallisto quant
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(prefix) - FASTQ Prefix
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   tuple
//     path(fa) - Reference FASTA
//     path(idx_files) - Index Files
//   val params - Additional Parameters
//   val out_dir - Output Directory
//   val shared_dir - Shared Output Directory
//
// output:
//   tuple => emit: quants
//     val(pat_name)
//     val(prefix)
//     val(dataset)
//     path('*abundance.tsv')
//   tuple => emit: h5_quants
//     val(pat_name)
//     val(prefix)
//     val(dataset)
//     path('*abundance.h5')
//   path('kallisto-*') - For publishDir

  label 'kallisto'
  publishDir "${out_dir}/${dataset}/${pat_name}/${prefix}", pattern: "kallisto-*", mode: 'copyNoFollow'
  publishDir "${shared_dir}/${dataset}/${pat_name}/${prefix}", pattern: "kallisto-*", mode: 'copyNoFollow'
  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  tuple path(fa), path(idx_files)
  val params
  val out_dir
  val shared_dir

  output:
  tuple val(pat_name), val(prefix), val(dataset), path('*abundance.tsv'), emit: quants
  tuple val(pat_name), val(prefix), val(dataset), path('*abundance.h5'), emit: h5_quants
  path('kallisto-*')

  script:
  """
  CMD="kallisto quant ${params} -t ${task.cpus} -i ${fa}.index -o . ${fq1} ${fq2}"

  VERSION=\$(kallisto | grep ^kallisto | cut -f 2 -d ' ')
  CMD_MD5=\$(echo \${CMD} | md5sum | cut -f 1 -d ' ' | cut -c1-5)
  FA_MD5=\$(md5sum ${fa} | cut -f 1 -d ' ' | cut -c1-5)
  FQ1_MD5=\$(md5sum ${fq1} | cut -f 1 -d ' ' | cut -c1-5)
  FQ2_MD5=\$(md5sum ${fq2} | cut -f 1 -d ' ' | cut -c1-5)

  OUTDIR="kallisto-\${VERSION}=\${CMD_MD5}/${fa}=\${FA_MD5}/${fq1}=\${FQ1_MD5}-${fq2}=\${FQ2_MD5}/"
  EXP_SHARED_DIR="${shared_dir}/${dataset}/${pat_name}/${prefix}/\${OUTDIR}/"


  mkdir -p \${OUTDIR}

  if [[ ! -d \${EXP_SHARED_DIR} ]]; then
    touch .real
    eval \${CMD}
    mv abundance.tsv ${dataset}-${pat_name}-${prefix}.abundance.tsv
    mv abundance.h5 ${dataset}-${pat_name}-${prefix}.abundance.h5
    ln -s \${PWD}/*abundance* \${OUTDIR}/
    ln -s .command.* \${OUTDIR}/
  else
    touch .copy
    cp -r \${EXP_SHARED_DIR}/* \${OUTDIR}
    cp -r \${EXP_SHARED_DIR}/* .
  fi
  """
}
