process mhcflurry_predict {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'mhcflurry_container'
  label 'mhcflurry'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/mhcflurry"
  cache 'lenient'

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(inp_file)
  path data_dir
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*mhcflurry.csv"), emit: mhcflurry_scores

  script:
  """
  mhcflurry-predict \
    ${inp_file} \
    --models ${data_dir} \
    ${parstr} \
    > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mhcflurry.csv
  """
}

process mhcflurry_make_input_file {

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'mhcflurry_container'
  label 'mhcflurry'
  cache 'lenient'

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(peptides), path(alleles)

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*mhcflurry_inp.tsv"), emit: mhcflurry_inpf

  script:
  """
  echo "peptide,allele" > ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mhcflurry_inp.tsv

  for i in `awk < ${peptides} '\$0 !~ /^>/'`;  do
    for j in `cat ${alleles} | sed 's/,/\\n/g'`;  do
     echo "\${i},\${j}" >> ${dataset}-${pat_name}-${norm_run}_${tumor_run}.mhcflurry_inp.tsv;
    done
  done
  """
}
