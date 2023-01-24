########################################
##### Process running on reference #####
########################################
process index_reference {
// Docstring
// Should explain inputs and output channels.

// require:                                       # Populates the process if loaded directly into main.nf.
//   params.<tool>$<input_file_or_parameter>      # File (e.g. reference FASTA. "params." prefix allows it to be captured by raft add-step command.
//   params.<tool>$<bwa_index_parameters          # String containing parameters passed by the user.  "params." prefix allows it to be captured by raft add-step command.

  storeDir "${params.shared_dir}/${fa}/bwa_index" # Directory where outputs are stored.
  tag "${fa}"                                     # Tag that identifies this process in Nextflow output.
  label '<tool>_container'                        # Nextflow label indicating which container to use.
  label '<tool>_<command>                            # Nextflow label indicating run parameters (e.g. CPU and RAM).

  input:                                          # Process inputs - Note that "require:" block and "input:" block have same number of elements.
  path fa                                         # Input file (e.g. FASTA). File will be referenced as ${fa} within "script:" block.
  val parstr                                      # Input string. parstr is standardized name for PARameter STRing within LENS.

  output                                          # Process outputs
  tuple path(fa), path("${fa}.*"), emit: idx_files# Tuple containing initially provided file and resulting index files. Note wildcard in "${fa}.*" -- This will capture all files that follow that pattern.

  script:                                         # "script:" block defining process behavior.
  """
  bwa index ${parstr} ${fa}                       # Literal command to be executed 
  """
}


##########################################################################
##### Process running on input(s) from single patient, single sample #####
##########################################################################
process sample_level_tool_single_run_name {
// Application: One patient, one sample (e.g. normal DNA ONLY, tumor RNA ONLY, etc.), single or multiple multiple input or output files
// The val(pat_name), val(run), val(dataset) pattern should be consistent when
// applicable (e.g. single patient, single sample, single or multiple input or
// output files).

// Docstring
// Should explain inputs and output channels.

// require:                                       # Populates the process if loaded directly into main.nf.
//   FQS                                          # User-provided parameter or channel. The CAPITAL LETTERS indiciate uers must change this. 
//   IDX_FILES                                    # User-provided parameter or channel. The CAPITAL LETTERS indiciate uers must change this. 
//   params.<tool>$<tool>_<command>_parameters       # String containing user-provided parameters. "params." prefix allows it to be captured by raft add-step command.


  publishDir "${params.output_dir}/${dataset}_${pat_name}_${run}/<tool>_<command> " # Directory where outputs are stored. Note the dataset -> patient -> run tagging of the output file.
  tag "${dataset}/${pat_name}/${run}"             # Tag that identifies this process in Nextflow output.
  label '<tool>_container'                        # Nextflow label indicating which container to use.
  label '<tool>_<command>                            # Nextflow label indicating run parameters (e.g. CPU and RAM).

  input:                                          # Process inputs - Note that "require:" block and "input:" block have same number of elements.
  tuple val(pat_name), val(run), val(dataset), path(fq1), path(fq2) # Input tuple defining sample-level inputs. Note that input files (fastqs) are "tagged" by the Patient_Name, Run_Name, and Dataset. It is crucial these "tagging" values maining this order in input and output strings.
  tuple path(fa), path(idx_files)                 # Input tuple defining reference files. 
  val parstr                                      # Input string. parstr is standardized name for PARameter STRing within LENS.

  output:                                         # Process outputs
  tuple val(pat_name), val(run), val(dataset), path('*.sam'), emit: sams # Tuple containing "tagged" `.sam` file. Note the "emit: sams" portion which allows this output so be called using <process_name>.out.sams by other processes/workflows.

  script:
  """
  <command>  ${parstr} ${fa} ${fq1} ${fq2} > ${dataset}-${pat_name}-${run}.sam -t ${task.cpus} # Literal command being run by process. Note that the output file name includes all tagging information (Patient_Name, Run_Name, and Dataset).
  """
}


#################################################################################
##### Process running on single input from single patient, multiple samples #####
#################################################################################
process sample_level_tool_multiple_run_names_single_input_file {
// Suitable for processes where you have a single input file per patient, but
// the file was generated using data from multiple samples from a patient (e.g.
// tumor and normal).
//
// The val(pat_name), val(norm_run), val(tumor_run), val(dataset) pattern
// should be consistent when applicable (e.g. single patient, multiple samples,
// one input or output).

// Docstring
// Should explain inputs and output channels.

// require: # Populates the process if loaded directly into main.nf.
//   VCFS  # User-provided parameter or channel. The CAPITAL LETTERS indiciate users must change this.
//   params.<tool>$<tool>_<command>_parameters       # String containing user-provided parameters. "params." prefix allows it to be captured by raft add-step command.

  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}" # Tag that identifies this process in Nextflow output.
  label '<tool>_container'                        # Nextflow label indicating which container to use.
  label '<tool>_<command>                            # Nextflow label indicating run parameters (e.g. CPU and RAM).
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/snpeff_ann" # Directory where outputs are stored. Note the dataset -> patient -> runs tagging of the output file.


  input:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path(vcf) # Input tuple defining sample-level inputs. Note that input files (fastqs) are "tagged" by the Patient_Name, Run_Name, and Dataset. It is crucial these "tagging" values maintain this order in input and output strings. Also, it is CRUCIAL the norm_run val occurs before the tumor_run within the tuple. 
  val parstr                                      # Input string. parstr is standardized name for PARameter STRing within LENS.

  output:
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*annot.vcf"), emit: annot_vcfs

  script:
  """
  /opt/snpeff/snpeff/bin/snpEff ann ${ref} ${vcf} -c ${ref}/snpEff.config > \${ORIG%.vcf*}.annot.vcf -dataDir \${PWD} # Literal command being executed
  """
}


####################################################################################
##### Process running on multiple inputs from single patient, multiple samples #####
####################################################################################
process sample_level_tool_multiple_run_names_multiple_input_files {
// Suitable for processes where patient's have multiple input files (one or
// more per sample). For example, a normal DNA BAM and BAI and  tumor DNA BAM
// and BAI used for variant calling.
//
// The val(pat_name), val(norm_run), val(tumor_run), val(dataset) pattern
// should be consistent when applicable (e.g. single patient, multiple samples,
// one input or output).

// require:  Populates the process if loaded directly into main.nf. 
//   NORM_TUMOR_BAMS_BAIS # User-provided parameter or channel. The CAPITAL LETTERS indiciate users must change this.
//   params.<tool>$<tool>_parameters # String containing user-provided parameters. "params." prefix allows it to be captured by raft add-step command.


  tag "${dataset}/${pat_name}/${norm_run}_${tumor_run}" # Tag that identifies this process in Nextflow output.
  label '<tool>_container'                        # Nextflow label indicating which container to use.
  label '<tool>_<command>                            # Nextflow label indicating run parameters (e.g. CPU and RAM).
  publishDir "${params.samps_out_dir}/${dataset}/${pat_name}/${norm_run}_${tumor_run}/snpeff_ann" # Directory where outputs are stored. Note the dataset -> patient -> runs tagging of the output file.


  input: # Process inputs
  tuple val(pat_name), val(dataset), val(norm_run), path(norm_bam), path(norm_bai), val(tumor_run), path(tumor_bam), path(tumor_bai) # Input tuple defining sample-level inputs. Note the normal sample is listed followed by normal file paths (bam and bai) followed by the tumor sample and its tumor file paths (bam and bai).
  val parstr # Input string. parstr is standardized name for PARameter STRing within LENS.

  output: # Process outputs 
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*snvs.vcf.gz"), emit: snv_vcfs # Tuple containing "tagged" `.snvs.vcf.gz` file. Note the "emit: snv_vcfs" portion which allows this output so be called using <process_name>.out.snv_vcfs by other processes/workflows.
  tuple val(pat_name), val(norm_run), val(tumor_run), val(dataset), path("*indels.vcf.gz"), emit: indel_vcfs # Tuple containing "tagged" `.indels.vcf.gz` file. Note the "emit: indel_vcfs" portion which allows this output so be called using <process_name>.out.indel_vcfs by other processes/workflows.


  script:
  """
  <command> ${parstr} --tumorBam ${tumor_bam} --normalBam ${norm_bam} --exome # Literal command being run by the process.
  """
}
