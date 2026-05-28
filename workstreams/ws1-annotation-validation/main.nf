#!/usr/bin/env nextflow

// WS1 skeleton: Input -> Conversion (T1) -> Core annotate+parse (T2) -> Exploration (T3).
// See workflow.mmd for the diagram this implements.

include { CONVERT  } from './modules/convert.nf'
include { ANNOTATE } from './modules/annotate.nf'
include { PARSE    } from './modules/parse.nf'
include { VALIDATE } from './modules/validate.nf'

workflow {
    if( !params.input )
        error "No input. Use --input <structure.pdb|.cif>  (or -profile test for the bundled 9CFN)"

    // Input
    structure_ch = Channel.fromPath(params.input, checkIfExists: true)

    // T1 Conversion: any structure -> standardized mmCIF
    mmcif_ch = CONVERT(structure_ch)

    // T2 Core: run each selected annotator on the mmCIF
    tools_ch = Channel.fromList( params.annotators.tokenize(',')*.trim() )
    raw_ch   = ANNOTATE( tools_ch.combine(mmcif_ch) )

    // T2 Parser: tool-native output -> standardized base-pairing mmCIF
    bp_ch = PARSE( raw_ch )

    // T3 Exploration: compare annotations across tools
    VALIDATE( bp_ch.map { tool, cif -> cif }.collect() )
}
