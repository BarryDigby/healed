#!/usr/bin/env nextflow

process kraken2_cmd {

  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  path db
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2), path(CMD), emit: cmd

  script:
  """
  CMD="kraken2 --db ${db} ${parstr} --paired ${fq1} ${fq2}"
  echo \${CMD} | tr -d '\n' > CMD
  """
}

process kraken2_uniq_dir {

  label 'kraken2_container'
  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2), path(CMD)
  path db

  output:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2), path(CMD), path(uniq_dir), path(meta), emit: total

  script:
  """
  # Actual CMD
  ACMD=\$(cat ${CMD})

  VERSION=\$(kraken2 --version | head -1 | cut -f 3 -d ' ')

  # Full checksums
  CMD_MD5=\$(echo \${ACMD} | md5sum | cut -f 1 -d ' ')
  DB_MD5=\$(md5sum ${db} | cut -f 1 -d ' ')
  FQ1_MD5=\$(md5sum ${fq1} | cut -f 1 -d ' ')
  FQ2_MD5=\$(md5sum ${fq2} | cut -f 1 -d ' ')

  # Truncated checksums
  CMD_MD5_T=\$(echo \${CMD_MD5} | cut -c1-5)
  DB_MD5_T=\$(echo \${CTAT_MD5} | cut -c1-5)
  FQ1_MD5_T=\$(echo \${FQ1_MD5} | cut -c1-5)
  FQ2_MD5_T=\$(echo \${FQ2_MD5} | cut -c1-5)

  OUTDIR="${dataset}/${pat_name}/${prefix}/kraken-\${VERSION}=\${CMD_MD5_T}/${fq1}=\${FQ1_MD5_T}-${fq2}=\${FQ2_MD5_T}/${db}=\${DB_MD5_T}"

  echo \${OUTDIR} | tr -d '\n' > uniq_dir

  echo "Dataset: ${dataset}" >> meta
  echo "Patient name: ${pat_name}" >> meta
  echo "Sample identifier: ${prefix}" >> meta
  echo "Container: ${task.container}" >> meta
  echo "Version: \${VERSION}" >> meta
  echo "Command: \${ACMD}" >> meta
  echo "Command checksum: \${CMD_MD5}" >> meta
  echo "Fastq1: ${fq1}" >> meta
  echo "Fastq1 checksum: \${FQ1_MD5}" >> meta
  echo "Fastq2: ${fq2}" >> meta
  echo "Fastq2 checksum: \${FQ2_MD5}" >> meta
  echo "Kraken database: ${db}" >> meta
  echo "Kraken database checksum: \${DB_MD5}" >> meta
  echo "Output directory: \${OUTDIR}" >> meta
  """
}

process kraken2_exec {

  label 'kraken2'
  tag "${dataset}/${pat_name}/${prefix}"
  publishDir "${proj_dir}/${uniq_dir.getText()}"
  storeDir "${shared_dir}/${uniq_dir.getText()}"

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2), path(CMD), val(uniq_dir), path(meta)
  path db
  val(proj_dir)
  val(shared_dir)

//  output:
//  tuple val(pat_name), val(prefix), val(dataset), path("*"), emit: kraken_outs

  script:
  """
  AUGMENTED_CMD="`cat ${CMD}` --threads ${task.cpus}"
  eval \${AUGMENTED_CMD}
  """
}

workflow kraken2 {
// require:
//   FQS
//   params.kraken2$kraken2_database
//   params.kraken2$kraken2_parameters
  take:
    fqs
    db
    parstr
  main:
    kraken2_cmd(
      fqs,
      db,
      parstr)
    kraken2_uniq_dir(
      kraken2_cmd.out.cmd,
      db)
    kraken2_exec(
      kraken2_uniq_dir.out.total,
      db,
      params.samps_out_dir,
      params.shared_dir)
//  emit:
//    kraken_outs = kraken2_exec.out.kraken_outs
}
