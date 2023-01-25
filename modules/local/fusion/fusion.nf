#!/usr/bin/env nextflow

include { trim_galore } from '../trim_galore/trim_galore.nf'
include { get_fastqs } from '../utilities/utilities.nf'
include { starfusion } from '../starfusion/starfusion.nf'
include { jstarfusion } from '../starfusion/starfusion.nf'

workflow manifest_to_starfusion_fusions {
// take:
//   ctat_ref - Trinity CTAT Reference
//   trim_galore_parameters - Trim Galore Parameters
//   starfusion_parameters - STARFusion Parameters
//   MANIFEST - RAFT Manifest
//
// emit:
//   fusions - STARFusion Fusions

// require:
//   params.fusion$manifest_to_starfusion_fusions$ctat_ref
//   params.fusion$manifest_to_starfusion_fusions$trim_galore_parameters
//   params.fusion$manifest_to_starfusion_fusions$starfusion_parameters
//   MANIFEST

  take:
    ctat_ref
    trim_galore_parameters
    starfusion_parameters
    manifest
  main:
    get_fastqs(
      manifest,
      params$fq_dir)
    raw_fqs_to_starfusion_fusions(
      ctat_ref,
      trim_galore_parameters,
      starfusion_parameters,
      get_fastqs.out.fastqs)
  emit:
    fusions = raw_fqs_to_starfusion_fusions.out.fusions
}


workflow raw_fqs_to_starfusion_fusions {
// take:
//   ctat_ref - Trinity CTAT Reference
//   trim_galore_parameters - Trim Galore Parameters
//   starfusion_parameters - STARFusion Parameters
//   FQS - Raw FASTQs
//
// emit:
//   fusions - STARFusion Fusions

// require:
//   params.fusion$manifest_to_starfusion_fusions$ctat_ref
//   params.fusion$raw_fqs_to_starfusion_fusions$trim_galore_parameters
//   params.fusion$raw_fqs_to_starfusion_fusions$starfusion_parameters
//   FQS

  take:
    ctat_ref
    trim_galore_parameters
    starfusion_parameters
    fqs
  main:
    trim_galore(
      fqs,
      trim_galore_parameters)
    procd_fqs_to_starfusion_fusions(
      ctat_ref,
      starfusion_parameters,
      trim_galore.out.procd_fqs)
  emit:
    fusions = procd_fqs_to_starfusion_fusions.out.fusions
}


workflow procd_fqs_to_starfusion_fusions {
// take:
//   ctat_ref - Trinity CTAT Reference
//   trim_galore_parameters - Trim Galore Parameters
//   starfusion_parameters - STARFusion Parameters
//   FQS - Processed FASTQs
//
// emit:
//   fusions - STARFusion Fusions

// require:
//   params.fusion$procd_fqs_to_starfusion_fusions$ctat_ref
//   params.fusion$procd_fqs_to_starfusion_fusions$starfusion_parameters
//   FQS

  take:
    ctat_ref
    starfusion_parameters
    fqs
  main:
    starfusion(
      fqs,
      ctat_ref,
      starfusion_parameters)
  emit:
    coding_effect_fusions = starfusion.out.coding_effect_fusions
    full_fusions = starfusion.out.full_fusions
}


workflow junctions_to_starfusion_fusions {
// take:
//   ctat_ref - Trinity CTAT Reference
//   starfusion_parameters - STARFusion Parameters
//   juncts - STAR Junctions
//
// emit:
//   coding_effect_fusions - STARFusion '--examine_coding_effect' Outputs
//   full_fusions - STARFusion Standard Fusions

// require:
//  params.fusion$procd_fqs_to_starfusion_fusions$ctat_ref
//  params.fusion$procd_fqs_to_starfusion_fusions$starfusion_parameters
//  JUNCTIONS

  take:
    ctat_ref
    starfusion_parameters
    juncts
  main:
    jstarfusion(
      juncts,
      ctat_ref,
      starfusion_parameters)
  emit:
    coding_effect_fusions = jstarfusion.out.coding_effect_fusions
    full_fusions = jstarfusion.out.full_fusions
}


process starfusion_fusions_to_extracted_reads {
// Extracts fusion-assocaited read identifiers from starfusion outputs.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path(fusion_data) - STARFusion fusions
//
// output:
//   tuple => emit: fusion_read_names
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) -  Dataset
//     path("*fusion_read_names") - Fusion-supporting Read Names

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fusion_data)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*fusion_read_names"), emit: fusion_read_names

  script:
  """
  JUNCTION_COL_IDX=`sed -n \$'1s/\\\\\t/\\\\\n/gp' ${fusion_data} | grep -nx 'JunctionReads' | cut -d: -f1`
  SPANNING_COL_IDX=`sed -n \$'1s/\\\\\t/\\\\\n/gp' ${fusion_data} | grep -nx 'SpanningFrags' | cut -d: -f1`

  #Grabbing neutral read identifier
  cut -f \$JUNCTION_COL_IDX,\$SPANNING_COL_IDX ${fusion_data} | \
  sed 's/\$/\\n/g' | \
  sed 's/,/\\n/g' | \
  sed 's/\\s/\\n/g' | \
  grep -v '^\$' | \
  cut -f 2 -d '@' > ${dataset}-${pat_name}-${run}.fusion_read_names

  #Grabbing first read identifier
  cut -f \$JUNCTION_COL_IDX,\$SPANNING_COL_IDX ${fusion_data} | \
  sed 's/\$/\\/1\\n/g' | \
  sed 's/,/\\/1\\n/g' | \
  sed 's/\\s/\\/1\\n/g' | \
  grep -v '^\$' | \
  cut -f 2 -d '@' >> ${dataset}-${pat_name}-${run}.fusion_read_names

  #Grabbing second read identifier
  cut -f \$JUNCTION_COL_IDX,\$SPANNING_COL_IDX ${fusion_data} | \
  sed 's/\$/\\/2\\n/g' | \
  sed 's/,/\\/2\\n/g' | \
  sed 's/\\s/\\/2\\n/g' | \
  grep -v '^\$' | \
  cut -f 2 -d '@' >> ${dataset}-${pat_name}-${run}.fusion_read_names
  """
}


process starfusion_fusions_to_fusion_txs {
// Extracts fusion-associated read identifiers from starfusion outputs.
//
// input:
//   tuple
//     val(pat_name) - Patient Name
//     val(run) -  Run Name
//     val(dataset) - Dataset
//     path(fusion_data) - STARFusion fusions
//
// output:
//   tuple => emit: fusion_transcripts
//     val(pat_name) -  Patient Name
//     val(run) -  Run Name
//     val(dataset) -  Dataset
//     path("*fusion_transcripts") - Transcripts with Fusions

  tag "${dataset}/${pat_name}/${run}"

  input:
  tuple val(pat_name), val(run), val(dataset), path(fusion_data)

  output:
  tuple val(pat_name), val(run), val(dataset), path("*fusion_transcripts"), emit: fusion_transcripts

  script:
  """
 cut -f 18,20 *tsv | \
  grep ENS | \
  sed 's/\\t/\\n/g' >> ${dataset}-${pat_name}-${run}.fusion_transcripts
  """
}
