include { sniffles } from '../sniffles/sniffles.nf'
include { cutesv } from '../cutesv/cutesv.nf'
include { svim_alignment } from '../svim/svim.nf'
include { delly_lr } from '../delly/delly.nf'

workflow bams_to_svs {
  take:
    bams
    aln_ref
    sv_caller_tool
    sv_caller_tool_parameters
  main:
    sv_caller_tool_parameters = Eval.me(sv_caller_tool_parameters)
    bams.view()
    if( sv_caller_tool =~ /sniffles/ ) {
      sniffles_parameters = sv_caller_tool_parameters['sniffles'] ? sv_caller_tool_parameters['sniffles'] : ''
      sniffles(
        bams,
        sniffles_parameters) 
    }
    if( sv_caller_tool =~ /cutesv/ ) {
      cutesv_parameters = sv_caller_tool_parameters['cutesv'] ? sv_caller_tool_parameters['cutesv'] : ''
      cutesv(
        bams,
        aln_ref,
        cutesv_parameters) 
    }
    if( sv_caller_tool =~ /svim/ ) {
      svim_parameters = sv_caller_tool_parameters['svim'] ? sv_caller_tool_parameters['svim'] : ''
      svim_alignment(
        bams,
        aln_ref,
        svim_parameters) 
    }
    if( sv_caller_tool =~ /delly/ ) {
      delly_parameters = sv_caller_tool_parameters['delly'] ? sv_caller_tool_parameters['delly'] : ''
      delly_lr(
        bams,
        aln_ref,
        delly_parameters) 
    }
//  emit:
}
