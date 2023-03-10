#!/usr/bin/env nextflow

/* Utilities.nf */
// Utilities for RAFT

process sideload_sample_files {
// Loads intermediate sample-level files. These files will typically serve as
// the starting point for a workflow (e.g. start with BAMs rather than FASTQs).
// This functionality is included for module testing purposes, but may be
// useful elsewhere as well.
// The FILE_REGEX variable should be a suffix associated with the desired file
// (e.g. ".bam", ".fastq.bz", ".vcf", etc.)

// require:
//   MANIFEST
//   FILE_REGEX
//   params.sideload_dir

  input:
  tuple val(pat_name), val(prefix), val(dataset)
  val suffix
  val sideload_dir

  output:
  tuple val(pat_name), val(prefix), val(dataset), path("${dataset}-${pat_name}-${prefix}${suffix}"), emit: sideloaded_files

  script:
  """
  discovered_file=\$(find -L ${sideload_dir} -name "${dataset}-${pat_name}-${prefix}${suffix}")
  ln -s \${discovered_file} \${PWD}/${dataset}-${pat_name}-${prefix}${suffix}
  """
}

process sideload_norm_tumor_files {
// Loads intermediate patient-level files derived from a normal/tumor pair of
// input files. These files will typically serve as the starting point for a
// workflow (e.g. start with BAMs rather than FASTQs).
// This functionality is included for module testing purposes, but may be
// useful elsewhere as well.
// The FILE_REGEX variable should be a suffix associated with the desired file
// (e.g. ".bam", ".fastq.bz", ".vcf", etc.)

// require:
//   MANIFEST
//   FILE_REGEX
//   params.sideload_dir

  input:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset)
  val suffix
  val sideload_dir

  output:
  tuple val(pat_name), val(norm_prefix), val(tumor_prefix), val(dataset), path("${dataset}-${pat_name}-${norm_prefix}*${tumor_prefix}${suffix}"), emit: sideloaded_files

  script:
  """
  discovered_file=\$(find -L ${sideload_dir} -name "${dataset}-${pat_name}-${norm_prefix}*${tumor_prefix}${suffix}")
  ln -s \${discovered_file} \${PWD}/${dataset}-${pat_name}-${norm_prefix}_${tumor_prefix}${suffix}
  """
}

process get_fastqs {
// Symlink and emit FASTQ pairs for a FASTQ Prefix of a Patient Name from a Dataset
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(prefix) - FASTQ Prefix
//     val(dataset) - Dataset
//     val(run_name) - Run Name
//     val fq_dir - Project-specific /fastq directory
//
// output:
//   tuple => emit: fastqs
//     val(pat_name) - Patient Name
//     val(run_name) - Sequence run name. Unique within dataset.
//     val(dataset) - Dataset
//     path("${prefix}*1*.f*q.gz") - FASTQ 1
//     path("${prefix}*2*.f*q.gz") - FASTQ 2

  tag "${dataset}/${pat_name}/${prefix}"

  input:
  tuple val(pat_name), val(run_name), val(dataset), val(prefix)
  val fq_dir

  output:
  tuple val(pat_name), val(run_name), val(dataset), path("*${run_name}_1*.f*q.gz"), path("*${run_name}_2*.f*q.gz"), emit: fastqs

  script:
  """
  fastq1=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_1\\|R1\\).*\\.\\(fq\\|fastq\\)\\.gz")
  fastq2=\$(find -L ${fq_dir} -regex ".*/${prefix}.*\\(_2\\|R2\\).*\\.\\(fq\\|fastq\\)\\.gz")
  ln -s \${fastq1} ${dataset}-${pat_name}-${run_name}_1.fastq.gz
  ln -s \${fastq2} ${dataset}-${pat_name}-${run_name}_2.fastq.gz
  """
}


process symlink_fastqs {
// Symlink FASTQs from RAFT global FASTQ directory into a project-specific
// FASTQ directory. Intended for usage where manifests are generated by dataset
// preparation scripts. This process emits an 'unlock' channel once symlinking
// as completted.
//
// input:
//  path manifest - Manifest generated by dataset preparation script
//  val global_fq_dir - RAFT Global FASTQ directory
//  val proj_fq_dir - Project-specific FASTQ directory
//  val filter - Filter for Sequencing_Method (if desired)

  input:
  path manifest
  val global_fq_dir
  val proj_fq_dir

  output:
  val 'true', emit: unlock

  script:
  """
#!/usr/bin/env python3

from glob import glob
import os
import re

dirs = []

def extract(f):
    prefixes = []
    with open(f) as ifo:
        hdr = ifo.readline()
        hdr = hdr.strip('\\n').strip().split('\t')
        prefix_col = hdr.index('File_Prefix')
        for row in ifo:
            row = row.strip('\\n').strip().split('\t')
            prefix = row[prefix_col]
            if prefix == 'NA':
                continue
            prefixes.append(prefix)
    return prefixes

def update_mounts_cfg(mounts_cfg, bound_dirs):
    print(mounts_cfg)
    print(bound_dirs)
    out = []
    with open(mounts_cfg, 'r') as ifo:
        line = ifo.readline()
        line = line.strip('\\n')
        paths = line.split(',')
        bind_dirs_to_add = []
        for bind_dir in bound_dirs:
            if not(any([bind_dir.startswith(path) for path in paths])):
                bind_dirs_to_add.append(bind_dir)
                for path in paths:
                    if path.startswith(bind_dir):
                        paths.remove(path)
        paths.extend(bind_dirs_to_add)
        paths = ','.join(paths) + '\\n'
        out.append(paths)

    with open(mounts_cfg, 'w') as fo:
        for row in out:
            fo.write(row)

prefixes = list(set(extract("${manifest}")))

for prefix in prefixes:
    #This is basically copying from RAFT functionality. Likely want to refactor
    #the RAFT script and have this as one of the functions we can call.
    globbed_dir = os.path.join("${global_fq_dir}", '**', prefix + '*')
    hits = glob(globbed_dir, recursive=True)
    for hit in hits:
        try:
            os.symlink(os.path.realpath(hit), os.path.join("${proj_fq_dir}", os.path.basename(hit)))
            dirs.append(os.path.dirname(os.path.realpath(hit)))
        except:
            pass

dirs = list(set(dirs))

update_mounts_cfg(os.path.join(os.path.dirname("${proj_fq_dir}"), 'workflow', 'mounts.config'), list(set(dirs)))
"""
}


process wait_signal_1 {
//  Process to hold for symlink_fastqs() to finish before emitting manifest.

  tag "${input}"

  input:
  val input
  val signal

  output:
  val input, emit: e_output

  shell:
  """
  echo "Got signal! Emitting output!"
  """
}


workflow extract_manifest_from_channel {
// Consume a manifest emitted by a channel and convert it to a format suitable
// for processing. This typically involves selecting only a handful of required
// columns and filtering based on Sequencing_Method.
//
// take:
//   raw_manifest -  Raw manifest emitted by a dataset preparation step.
//   molecule_filter - A regular expression used for filtering samples (e.g. '^RNA')
//   separator - Delimiter used to parse rows of raw_manifest
//
// emit:
//   manifest - Formatted manifest ready for processing (by manifest_to_* steps)

// require:
//   MANIFEST
//   params.utilities$extract_manifest_from_channel$molecule_filter
//   params.utilities$extract_manifest_from_channel$separator
  take:
    raw_manifest
    molecule_filter
    separator
  main:
    symlink_fastqs(raw_manifest, params.global_fq_dir, params.fq_dir)

    raw_manifest.splitCsv(header: true, sep: separator)
    .map{ row -> tuple("${row.Patient_Name}", "${row.Run_Name}", "${row.Dataset}", "${row.File_Prefix}", "${row.Sequencing_Method}", "${row.Normal}") }
    .filter{ it[4] =~ /${molecule_filter}/ }
    .set{ manifest }

    wait_signal_1(manifest.toList(), symlink_fastqs.out.unlock)
  emit:
    manifest = wait_signal_1.out.e_output.flatMap{ it -> it }
}


workflow filter_manifest_from_channel {
// TODO: Change variable to be manifest_channel.
// Consume a manifest emitted by a channel and convert it to a format suitable
// for processing. This typically involves selecting only a handful of required
// columns and filtering based on Sequencing_Method.
//
// take:
//   raw_manifest -  Raw manifest emitted by a dataset preparation step.
//   filter - A regular expression used for filtering samples (e.g. '^RNA')
//   separator - Delimiter used to parse rows of raw_manifest
//
// emit:
//   manifest - Formatted manifest ready for processing (by manifest_to_* steps)

// require:
//  params.utilities$filter_manifest_from_channel$manifest
//  params.utilities$filter_manifest_from_channel$filter
//  params.utilities$filter_manifest_from_channel$separator
  take:
    manifest
    filter
    separator
  main:
    manifest
    .filter{ it[3] =~ /${filter}/ }
    .set{ filtered_manifest }

  emit:
    filtered_manifest = filtered_manifest
}


workflow extract_manifest_from_file {
// TODO: Change variable to be manifest_file.
// Consume a manifest from a path and convert it to a format suitable for
// processing. This typically involves selecting only a handful of required
// columns and filtering based on Sequencing_Method. This workflow should
// commonly be called on the project-wide manifest
// (/path/to/raft/projects/<project>/metadata/<project>_manifest.csv)
//
// take:
//   raw_manifest - Path containing raw manifest.
//   molecule_filter - A regular expression used for filtering samples (e.g. '^RNA')
//   separator - Delimiter used to parse rows of raw_manifest
//
// emit:
//   manifest - Formatted manifest ready for processing (by manifest_to_* steps)

// require:
//   MANIFEST
//   params.utilities$extract_manifest_from_file$molecule_filter
//   params.utilities$extract_manifest_from_file$separator
  take:
    raw_manifest
    molecule_filter
    separator
  main:
    symlink_fastqs(raw_manifest, params.global_fq_dir, params.fq_dir)

    raw_manifest.splitCsv(header: true, sep: separator)
    .map{ row -> tuple("${row.Patient_Name}", "${row.Run_Name}", "${row.Dataset}", "${row.File_Prefix}", "${row.Sequencing_Method}", "${row.Normal}") }
    .filter{ it[4] =~ /${molecule_filter}/ }
    .set{ manifest }

    wait_signal_1(manifest.toList(), symlink_fastqs.out.unlock)
  emit:
    manifest = wait_signal_1.out.e_output.flatMap{ it -> it }
}


workflow combine_sample_files {
// Combines multiple files derived from 
// require:
//  SET_1
//  SET_2
  take:
    set_1
    set_2
  main:
    set_1.combine(set_2, by: [0, 1, 2]).set{ combined_set }
  emit:
    combined_set
}


workflow combine_sample_somatic_files {
// require:
//  SET_1
//  SET_2
  take:
    set_1
    set_2
  main:
    set_1.combine(set_2, by: [0, 1, 2, 3]).set{ combined_set }
  emit:
    combined_set
}

workflow combine_patient_samples {
// require:
//  TOTAL_SET
//  FILTER_1
//  FILTER_2
//  MANIFEST
  take:
    total_set
    filter_1
    filter_2
    manifest
  main:
    manifest.filter{ it[5] =~ filter_1 }.set{ set_1 }
    manifest.filter{ it[5] =~ filter_2 }.set{ set_2 }
    total_set.join(set_1, by: [0, 1, 2]).map{ [it[0], it[1], it[2], it[3], it[4]] }.set{ actual_set_1 }
    total_set.join(set_2, by: [0, 1, 2]).map{ [it[0], it[1], it[2], it[3], it[4]] }.set{ actual_set_2 }
    actual_set_1.combine(actual_set_2, by: [0, 2]).set{ combined_set }
  emit:
    combined_set
}

workflow combine_rel_patient_samples {
// This  is for dealing with samples that have gone through some sort of
// relative analysis that still produces normal and tumor outputs (e.g. ABRA2).
// require:
//  TOTAL_SET
  take:
    total_set
  main:
    total_set.filter{ it[1] =~ "^n" }.set{ norm_set }
    total_set.filter{ it[1] =~ "^a" }.set{ tumor_set }
    norm_set.map{ [it[0], it[1].split('-rel-')[1], it[1].split('-rel-')[0], it[2], it[3], it[4]] }.set{ remapped_norm_set }
    tumor_set.map{ [it[0], it[1].split('-rel-')[1], it[1].split('-rel-')[0], it[2], it[3], it[4]] }.set{ remapped_tumor_set }
    remapped_norm_set.join(remapped_tumor_set, by: [0, 1, 3]).map{ [it[0], it[2], it[3], it[4], it[5], it[6], it[7], it[8]] }.set{ combined_rel_set }
  emit:
    combined_rel_set
}

workflow make_norm_tumor_prefixes {
// require:
//  FILTER_1
//  FILTER_2
//  MANIFEST
  take:
    filter_1
    filter_2
    manifest
  main:
    manifest.filter{ it[4] =~ filter_1 }.set{ set_1 }
    manifest.filter{ it[4] =~ filter_2 }.set{ set_2 }
    set_1.join(set_2, by: [0, 2]).map{ [it[0], it[2], it[5], it[1]] }.set{ norms_tumors }
  emit:
    norms_tumors
}

workflow filter_channel_by_manifest {
// require:
//   INIT_CHANNEL
//   FILTER
//   FILTER_IDX
//   MANIFEST
  take:
    init_channel
    filt_str
    filt_idx
    manifest
  main:
    manifest.filter{ it[filt_idx] =~ filt_str}.map{ [it[0], it[1], it[2]] }.set{ filt_manifest }
    init_channel.join(filt_manifest, by: [0, 1, 2]).set{ filt_channel }
  emit: 
    filt_channel
} 
