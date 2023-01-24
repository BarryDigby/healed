#!/usr/bin/env nextflow

process netmhcpan {
// Predicts binding of peptides to any MHC molecule of known sequence using
// artificial neural networks (ANNs).
//
// input:
// output:

// require:
//   PEPTIDES_AND_ALLELES

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'netmhcpan_container'
  label 'netmhcpan'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/netmhcpan"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(peptides), path(alleles)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*netmhcpan*"), emit: c1_antis

  script:
  """
  export TMPDIR="\${PWD}"
  /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan -listMHC > allowed_alleles

  for i in `cat *netmhcpan.alleles* | sed 's/,/\\n/g'`; do echo \${i}; grep "^\${i}\$" allowed_alleles | cat >> final_alleles; done

  EXPLICIT_ALLELES=\$(cat final_alleles | sed -z 's/\\n/,/g ' | sed 's/,\$//g')
  PEPTIDES_FILENAME="${peptides}"
  
  if [[ -s ${peptides} ]]; then
    /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan -BA -s -a \${EXPLICIT_ALLELES} -inptype 0 -f ${peptides} > \${PEPTIDES_FILENAME%.aa.fa}.netmhcpan.txt
  else
    echo "NO PEPTIDES" > \${PEPTIDES_FILENAME%.aa.fa}.netmhcpan.txt
  fi
  """
}



process netmhcpan_rna {
// Predicts binding of peptides to any MHC molecule of known sequence using
// artificial neural networks (ANNs).
//
// input:
// output:

// require:
//   PEPTIDES_AND_ALLELES

  tag "${dataset}/${pat_name}/${run}"
  label 'netmhcpan_container'
  label 'netmhcpan'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/netmhcpan"

  input:
  tuple val(pat_name), val(run), val(dataset), path(peptides), path(alleles)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*netmhcpan*"), emit: c1_antis

  script:
  """
  export TMPDIR="\${PWD}"
  /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan -listMHC > allowed_alleles

  for i in `cat *netmhcpan.alleles* | sed 's/,/\\n/g'`; do echo \${i}; grep "^\${i}\$" allowed_alleles | cat >> final_alleles; done

  EXPLICIT_ALLELES=\$(cat final_alleles | sed -z 's/\\n/,/g ' | sed 's/,\$//g')
  PEPTIDES_FILENAME="${peptides}"
  if [[ -s ${peptides} ]]; then
    export TMPDIR="\${PWD}"
    /netMHCpan-4.1/Linux_x86_64/bin/netMHCpan -BA -s -a \${EXPLICIT_ALLELES} -inptype 0 -f ${peptides} > \${PEPTIDES_FILENAME%.aa.fa}.netmhcpan.txt
  else
    PEPTIDES_FILENAME="${peptides}"
    echo "NO PEPTIDES" > \${PEPTIDES_FILENAME%.aa.fa}.netmhcpan.txt
  fi
  """
}
