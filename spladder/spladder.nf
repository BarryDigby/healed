process spladder_build{

  tag "${dataset}/${pat_name}/${run}"
  label 'spladder_container'
  label 'spladder_build'

  input:
  tuple val(pat_name), val(run), val(dataset), path(bam)
  path gtf
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*_graph.json"), emit: graph_jsons

  script:
  """
  spladder build \
    --gtf ${gtf} \
    --bams ${bam} \
    ${parstr}
  """
}
