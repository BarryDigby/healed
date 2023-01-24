#!/usr/bin/env nextflow

process netmhcstabpan {
// Predicts binding stability of peptides to any known MHC molecule using
// artificial neural networks (ANNs)
//
// input:
// output:

// require:
//   PEPTIDES_AND_ALLELES

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}"
  label 'netmhcstabpan_container'
  label 'netmhcstabpan'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/netmhcstabpan"

  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(peptides), path(alleles)
  val parstr

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*netmhcstabpan*"), emit: c1_stabs

  script:
  """
  export TMPDIR="\${PWD}"
  /netMHCstabpan-1.0/bin/netMHCstabpan -listMHC > allowed_alleles                             
                                                                                                    
  for i in `cat *netmhcpan.alleles* | sed 's/,/\\n/g'`; do echo \${i}; grep "^\${i}\$" allowed_alleles | cat >> final_alleles; done

  EXPLICIT_ALLELES=\$(cat final_alleles | sed -z 's/\\n/,/g ' | sed 's/,\$//g')                                                                                                  
  PEPTIDES_FILENAME="${peptides}"
  
  if [[ -s ${peptides} ]]; then
    /netMHCstabpan-1.0/bin/netMHCstabpan -l 8,9,10,11 -a \${EXPLICIT_ALLELES} -inptype 0 -f ${peptides} > \${PEPTIDES_FILENAME%.aa.fa}.netmhcstabpan.txt
  else
    echo "NO PEPTIDES" > \${PEPTIDES_FILENAME%.aa.fa}.netmhcstabpan.txt
  fi
  """
}


process netmhcstabpan_rna {
// Predicts binding stability of peptides to any known MHC molecule using
// artificial neural networks (ANNs)
//
// input:
// output:

// require:
//   PEPTIDES_AND_ALLELES

  tag "${dataset}/${pat_name}/${run}"
  label 'netmhcstabpan_container'
  label 'netmhcstabpan'
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${run}/netmhcstabpan"

  input:
  tuple val(pat_name), val(run), val(dataset), path(peptides), path(alleles)
  val parstr

  output:
  tuple val(pat_name), val(run), val(dataset), path("*netmhcstabpan*"), emit: c1_stabs

  script:
  """
  export TMPDIR="\${PWD}"
  /netMHCstabpan-1.0/bin/netMHCstabpan -listMHC > allowed_alleles                             
                                                                                                    
  for i in `cat *netmhcpan.alleles* | sed 's/,/\\n/g'`; do echo \${i}; grep "^\${i}\$" allowed_alleles | cat >> final_alleles; done
                                                                                                    
  EXPLICIT_ALLELES=\$(cat final_alleles | sed -z 's/\\n/,/g ' | sed 's/,\$//g') 
  PEPTIDES_FILENAME="${peptides}"
  if [[ -s ${peptides} ]]; then
    /netMHCstabpan-1.0/bin/netMHCstabpan -l 8,9,10,11 -a \${EXPLICIT_ALLELES} -inptype 0 -f ${peptides} > \${PEPTIDES_FILENAME%.aa.fa}.netmhcstabpan.txt
  else
    PEPTIDES_FILENAME="${peptides}"
    echo "NO PEPTIDES" > \${PEPTIDES_FILENAME%.aa.fa}.netmhcstabpan.txt
  fi
  """
}
