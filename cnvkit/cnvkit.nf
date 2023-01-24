process cnvkit_batch {

  input:
    path bams
    path bais
    path targets
    path fa
  output:
    path "cnvkit.outputs", emit: cnvkit_outputs
  script:
  """
  cnvkit.py batch *-ad*bam --normal *-nd*bam \
      --targets ${targets} \
      --fasta ${fa} \
      --output-reference joint_normals.cnn \
      --output-dir cnvkit.outputs
  """
}
