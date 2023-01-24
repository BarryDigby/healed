process tximport {
  
  label "tximport_container"
  label "tximport"

  input:
  tuple val(pat_name), val(run), val(dataset), path(quants)
  path tx2gene
  val quant_type
  val suffix

  output:
  tuple val(pat_name), val(run), val(dataset), path("*txi*${suffix}"), emit: gene_count_mtxs

  script:
  """
  echo "library(tximport)" > CMD
  echo "tx2gene <- read.table(\\${params.tx2gene}\\", header=T, sep=',')" >> CMD
  echo "quants <- list.files(pattern=\\"*sf\\", all.files=T)" >> CMD
  echo "names(quants) <- quants" >> CMD
  echo "txi <- tximport(files=quants, type=\\"${params.quant_type}\\", tx2gene=tx2gene)" >> CMD
  echo "write.table(txi, \\"${dataset}-${pat_name}-${run}txi${suffix}\\", col.names=T, row.names=T, sep="\t", quote=F)" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}
