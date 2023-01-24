#!/usr/bin/env nextflow

process antigen_garnish_foreignness {
// Runs antigen.garnish foreignness calculation
//
// input:
//
// output:

// require:

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'antigen_garnish_container'
  label 'antigen_garnish_foreignness'


  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(peptides_file)
  path ag_data_dir
  val species

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*ag_foreign.tsv"), emit: foreignness_files

  script:
  """
  export PATH=\$PATH:"${ag_data_dir}"
  export AG_DATA_DIR="${ag_data_dir}"
  chmod 770 antigen.garnish/blastp
  chmod +x antigen.garnish/blastp
  echo "#!/usr/bin/env Rscript" > CMD
  echo "library(magrittr)" >> CMD
  echo "library(data.table)" >> CMD
  echo "library(antigen.garnish)" >> CMD
  echo "peptide_table <- read.table(file = '${peptides_file}', header=TRUE, sep='\\t')" >> CMD
  echo "peptides <- as.vector(peptide_table\\\$Peptide)" >> CMD
  echo "foreigns <- peptides %>% foreignness_score(db = '${species}')" >> CMD
  echo "write.table(foreigns, file = '${dataset}-${pat_name}-${norm_run}_${tumor_run}.ag_foreign.tsv', row.names=FALSE, quote=FALSE, sep=',')" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}


process antigen_garnish_foreignness_rna {
// Runs antigen.garnish foreignness calculation
//
// input:
//
// output:

// require:

  tag "${dataset}/${pat_name}/${run}"
  label 'antigen_garnish_container'
  label 'antigen_garnish_foreignness'


  input:
  tuple val(pat_name), val(run), val(dataset), path(peptides_file)
  path ag_data_dir
  val species

  output:
  tuple val(pat_name), val(run), val(dataset), path("*ag_foreign.tsv"), emit: foreignness_files

  script:
  """
  export PATH=\$PATH:"${ag_data_dir}"
  export AG_DATA_DIR="${ag_data_dir}"
  chmod 770 antigen.garnish/blastp
  chmod +x antigen.garnish/blastp
  echo "#!/usr/bin/env Rscript" > CMD
  echo "library(magrittr)" >> CMD
  echo "library(data.table)" >> CMD
  echo "library(antigen.garnish)" >> CMD
  echo "peptide_table <- read.table(file = '${peptides_file}', header=TRUE, sep='\\t')" >> CMD
  echo "peptides <- as.vector(peptide_table\\\$Peptide)" >> CMD
  echo "foreigns <- peptides %>% foreignness_score(db = '${species}')" >> CMD
  echo "write.table(foreigns, file = '${dataset}-${pat_name}-${run}.ag_foreign.tsv', row.names=FALSE, quote=FALSE, sep=',')" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}

process antigen_garnish_dissimilarity {
// Runs antigen.garnish dissimilarity calculation
//
// input:
//
// output:

// require:

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'antigen_garnish_container'
  label 'antigen_garnish_dissimilarity'


  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path(peptides_file)
  path ag_data_dir
  val species

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(rna_run), val(dataset), path("*ag_dissim.tsv"), emit: dissimilarity_files

  script:
  """
  export PATH=\$PATH:"${ag_data_dir}"
  export AG_DATA_DIR="${ag_data_dir}"
  chmod 770 antigen.garnish/blastp
  chmod +x antigen.garnish/blastp
  echo "#!/usr/bin/env Rscript" > CMD
  echo "library(magrittr)" >> CMD
  echo "library(data.table)" >> CMD
  echo "library(antigen.garnish)" >> CMD
  echo "peptide_table <- read.table(file = '${peptides_file}', header=TRUE, sep='\\t')" >> CMD
  echo "peptides <- as.vector(peptide_table\\\$Peptide)" >> CMD
  echo "dissims <- peptides %>% dissimilarity_score(db = '${species}')" >> CMD
  echo "write.table(dissims, file = '${dataset}-${pat_name}-${norm_run}_${tumor_run}.ag_dissim.tsv', row.names=FALSE, quote=FALSE, sep=',')" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}


process antigen_garnish_dissimilarity_rna {
// Runs antigen.garnish dissimilarity calculation
//
// input:
//
// output:

// require:

  tag "${dataset}/${pat_name}/${run}"
  label 'antigen_garnish_container'
  label 'antigen_garnish_dissimliarity'


  input:
  tuple val(pat_name), val(run), val(dataset), path(peptides_file)
  path ag_data_dir
  val species

  output:
  tuple val(pat_name), val(run), val(dataset), path("*ag_dissim.tsv"), emit: dissimilarity_files

  script:
  """
  export PATH=\$PATH:"${ag_data_dir}"
  export AG_DATA_DIR="${ag_data_dir}"
  chmod 770 antigen.garnish/blastp
  chmod +x antigen.garnish/blastp
  echo "#!/usr/bin/env Rscript" > CMD
  echo "library(magrittr)" >> CMD
  echo "library(data.table)" >> CMD
  echo "library(antigen.garnish)" >> CMD
  echo "peptide_table <- read.table(file = '${peptides_file}', header=TRUE, sep='\\t')" >> CMD
  echo "peptides <- as.vector(peptide_table\\\$Peptide)" >> CMD
  echo "dissims <- peptides %>% dissimilarity_score(db = '${species}')" >> CMD
  echo "write.table(dissims, file = '${dataset}-${pat_name}-${run}.ag_dissim.tsv', row.names=FALSE, quote=FALSE, sep=',')" >> CMD
  chmod +x CMD
  Rscript CMD
  """
}
