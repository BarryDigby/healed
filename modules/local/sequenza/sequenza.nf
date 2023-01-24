#!/usr/bin/env nextflow

include { extract_chroms_from_bed } from '../ref_utils/ref_utils.nf'

process sequenza_gc_wiggle {
// require:
//   REFERENCE
//   params.sequenza$sequenza_gc_wiggle_window_size
//   params.sequenza$sequenza_gc_wiggle_parameters

  tag "${fa}"
  label 'sequenza_container'
  label 'sequenza_wiggle'

  input:
  path fa
  val window
  val parstr

  output:
  tuple path(fa), path("*wig.gz"), emit: gc_wig

  script:
  """
  sequenza-utils gc_wiggle -w ${window} --fasta ${fa} -o gc${window}Base.wig.gz
  """
}



process sequenza_bam2seqz_sub {
// require:
//   NORM_TUMOR_BAMS_BAIS
//   FA_W_GC

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'sequenza_container'
  label 'sequenza_bam2seqz'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/sequenza_bam2seqz"

  input:
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai)
  tuple path(fa), path(gc_wig_gz)
  path chroms_list

  output:
  // Cleaning up order of output tuple...
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.seqz.gz"), emit: seqz

  script:
  """
  RAW_CHROMS=\$(cat ${chroms_list} | grep -v chrM | grep -v chrY)
  CHROMS="\${RAW_CHROMS}"
  sequenza-utils bam2seqz -n ${norm_bam} -t ${tumor_bam} --chromosome \${CHROMS} --fasta ${fa} -gc ${gc_wig_gz} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.seqz.gz --parallel ${task.cpus}
  """
}


workflow sequenza_bam2seqz {
 take:
   bams_bais
   fa_w_gc
   bed
  main:
    extract_chroms_from_bed(
      bed)
    sequenza_bam2seqz_sub(
      bams_bais,
      fa_w_gc,
      extract_chroms_from_bed.out.chroms_list)
    merge_seqz(
      sequenza_bam2seqz_sub.out.seqz)
  emit:
    seqz = merge_seqz.out.merged_seqz
}


process sequenza_seqz_binning {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'sequenza_container'
  label 'sequenza_seqz_binning'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/sequenza_seqz_binning"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(seqz)
  val window

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*.small.seqz.gz"), emit: small_seqz

  script:
  """
  sequenza-utils seqz_binning --seqz ${seqz} -w ${window} -o ${dataset}-${pat_name}-${norm_run}_${tumor_run}.small.seqz.gz
  """
}



process sequenza_extract {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'sequenza_container'
  label 'sequenza_extract'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/sequenza_extract"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(small_seqz)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*extracted.Rdata"), emit: extracts

  script:
  """
  echo #!/usr/bin/env Rscript" >> CMD
  echo "library(sequenza)" >> CMD
  echo "extracted <- sequenza.extract(\\"${small_seqz}\\")" >> CMD
  echo "saveRDS(extracted, file=\\"${dataset}-${pat_name}-${norm_run}_${tumor_run}.extracted.Rdata\\")" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}


process sequenza_fit {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'sequenza_container'
  label 'sequenza_fit'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/sequenza_fit"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(extracted_rdata)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(extracted_rdata), path("*fit.Rdata"), emit: fits

  script:
  """
  echo "#!/usr/bin/env Rscript" >> CMD
  echo "library(sequenza)" >> CMD
  echo "extracted = readRDS(\\"${extracted_rdata}\\")" >> CMD
  echo "CP <- sequenza.fit(extracted)" >> CMD
  echo "saveRDS(CP, file=\\"${dataset}-${pat_name}-${norm_run}_${tumor_run}.fit.Rdata\\")" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}


process sequenza_result {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'sequenza_container'
  label 'sequenza_result'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/sequenza_result"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(extracted), path(fit)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*sequenza.results"), emit: results
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*segment*"), emit: segments
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*alternative_solutions*"), emit: solutions
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*segment*"), path("*alternative_solutions*"), emit: segments_and_solutions

  script:
  """
  echo "#!/usr/bin/env Rscript" >> CMD
  echo "library(sequenza)" >> CMD
  echo "extracted =  readRDS(\\"${extracted}\\")" >> CMD
  echo "fit = readRDS(\\"${fit}\\")" >> CMD
  #Keeping the 1:19 line here. This is required for Sequenza to successfully run on mouse samples. Ideally, Sequenza should be provided the chromosome list from something like samtools view of the input bams.
  #echo "sequenza.results(sequenza.extract = extracted, cp.table = fit, sample.id = \\"${dataset}-${pat_name}-${norm_run}_${tumor_run}\\", out.dir=\\"${dataset}-${pat_name}-${norm_run}_${tumor_run}.sequenza.results\\", chromosome.list=1:19)" >> CMD
  echo "sequenza.results(sequenza.extract = extracted, cp.table = fit, sample.id = \\"${dataset}-${pat_name}-${norm_run}_${tumor_run}\\", out.dir=\\"${dataset}-${pat_name}-${norm_run}_${tumor_run}.sequenza.results\\")" >> CMD
  chmod +x CMD
  Rscript CMD
  cp *sequenza.results/*segment* .
  cp *sequenza.results/*solution* .
  """
}


process merge_seqz {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  container 'starkngs/bgzip:latest'
//  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/merge_seqz"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(seqz) 

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*merged.seqz.gz"), emit: merged_seqz

  script:
  """
  SEQZ=\$(ls -v *seqz.gz)
  zcat \${SEQZ} | gawk '{if (NR!=1 && \$1 != "chromosome") {print \$0}}' | bgzip > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.merged.seqz.gz
  """
}
