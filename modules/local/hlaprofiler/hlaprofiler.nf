#!/usr/bin/env nextflow

process hlaprofiler_predict {
// Runs HLAProfiler.pl predict
// Uses the database included with the Docker image.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path(fq1) - FASTQ 1
//     path(fq2) - FASTQ 2
//   val parstr - Additional Parameters
//
// output:
//   tuple => emit: calls
//     val(pat_name) - Patient Name
//     val(run) - Run Name
//     val(dataset) - Dataset
//     path("*HLATypes.txt") - Calls file
//
// require:
//   FQS
//   params.hlaprofiler$hlaprofiler_predict$parameters

  tag "${dataset}/${pat_name}/${run}"
  label 'hlaprofiler_container'
  label 'hlaprofiler_predict'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/hlaprofiler_predict"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*HLATypes.txt"), emit: calls

  script:
  """
  perl /HLAProfiler-v.1.10beta/bin/HLAProfiler.pl predict \
    -fastq1 ${fq1} \
    -fastq2 ${fq2} \
    -database_name hla_database \
    -database_dir /opt/HLAProfiler-1.0.0-db_only \
    -reference /opt/HLAProfiler-1.0.0-db_only/hla_database/data/reference/hla.ref.merged.fa \
    -kraken_path /HLAProfiler-v.1.10beta/kraken-0.10.5-beta-ea.1 \
    -output_dir subwork \
    -allele_refinement all \
    -l sample.HLAProfiler.log \
    ${parstr} \
    -threads ${task.cpus}
  mv subwork/*/*HLATypes.txt ${dataset}-${pat_name}-${run}.HLATypes.txt
  find subwork -name "*fq" -exec gzip {} \\;
  """
}
