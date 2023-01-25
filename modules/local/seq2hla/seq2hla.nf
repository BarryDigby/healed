process seq2hla {

  tag "${dataset}/${pat_name}/${prefix}"

  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${prefix}/seq2hla"
  label 'seq2hla_container'
  label 'seq2hla'
  cache 'lenient'

  input:
  tuple val(pat_name), val(prefix), val(dataset), path(fq1), path(fq2)
  val parstr

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("*ClassI.HLAgenotype4*"), emit: calls

  script:
  """
  python /download/seq2hla/seq2HLA.py \
    -1 ${fq1} \
    -2 ${fq2} \
    -r ${dataset}-${pat_name}-${prefix} \
    -p ${task.cpus} \
    ${parstr}
  """
}
