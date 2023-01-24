#!/usr/bin/env nextflow

process generic_process_1_out {
  
  input:
    path inputs
    val process_container
    val num_cpus
    val amt_ram
    val outdir
    val out_regex_1
    val cmd
    
  container "${process_container}"
  publishDir "${outdir}"
  cpus "${num_cpus}"
  memory "${amt_ram}"
  
  output:
    path out_regex_1, emit: path_1

  script:
  """
  ${cmd}
  """
}


process generic_process_2_out {
  
  input:
    path inputs
    val process_container
    val num_cpus
    val amt_ram
    val outdir
    val out_regex_1
    val out_regex_2
    val cmd
    
  container "${process_container}"
  publishDir "${outdir}"
  cpus "${num_cpus}"
  memory "${amt_ram}"
  
  output:
    path out_regex_1, emit: path_1
    path out_regex_2, emit: path_2

  script:
  """
  ${cmd}
  """
}



process generic_process_3_out {
  
  input:
    path inputs
    val process_container
    val num_cpus
    val amt_ram
    val outdir
    val out_regex_1
    val out_regex_2
    val out_regex_3
    val cmd
  
  container "${process_container}"
  publishDir "${outdir}"
  cpus "${num_cpus}"
  memory "${amt_ram}"
  
  output:
    path out_regex_1, emit: path_1
    path out_regex_2, emit: path_2
    path out_regex_3, emit: path_3
  
  script:
  """
  ${cmd}
  """
}


process generic_process_4_out {
  
  input:
    path inputs
    val process_container
    val num_cpus
    val amt_ram
    val outdir
    val out_regex_1
    val out_regex_2
    val out_regex_3
    val out_regex_4
    val cmd
    
  container "${process_container}"
  publishDir "${outdir}"
  cpus "${num_cpus}"
  memory "${amt_ram}"

  output:
    path out_regex_1, emit: path_1
    path out_regex_2, emit: path_2
    path out_regex_3, emit: path_3
    path out_regex_4, emit: path_4
  
  script:
  """
  ${cmd}
  """
}

process tagged_generic_process {

  container "${process_container}"
  publishDir "${outdir}"

  input:
  tuple val(samp_id), path(inputs)
  val process_container
  val outdir
  val emit_id
  val out_regex
  val cmd

  output:
  tuple val(samp_id), path("${out_regex}"), emit: emit_id

  script:
  """
  ${cmd}
  """
}
