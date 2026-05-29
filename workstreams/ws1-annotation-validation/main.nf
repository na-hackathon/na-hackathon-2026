#!/usr/bin/env nextflow

// End-to-end annotation pipeline:
//   structure (.pdb|.cif) -> T1 CONVERT (MaxIT) -> T2 ANNOTATE (fr3d|rnapolis)
//   -> T2 PARSE -> standardized base-pairing mmCIF (per tool)
//   -> T3 VALIDATE (cross-tool comparison) + optional VISUALIZE (--visualize: per-tool graphs)
//
// Usage:
//   nextflow run main.nf <structure.pdb|.cif>                              # both tools, default
//   nextflow run main.nf --input <file> --annotators fr3d                  # single tool
//   nextflow run main.nf --input <file> --annotators fr3d,rnapolis         # both tools + validate
//
// The per-tool annotated mmCIF is published to ${params.outdir}/basepairs/${tool}_basepairs.cif.

include { CONVERT  } from './modules/convert.nf'
include { ANNOTATE } from './modules/annotate.nf'
include { PARSE    } from './modules/parse.nf'
include { VALIDATE; VISUALIZE } from './modules/validate.nf'

workflow {
    def input = params.input ?: ( args ? args[0] : null )
    if( !input )
        error "No input. Pass a structure positionally (nextflow run main.nf <file>) or with --input <file>  (or -profile test for the bundled 9CFN)"

    def tools = params.annotators.tokenize(',')*.trim().findAll { it }
    if( !tools )
        error "No annotators selected. Set --annotators fr3d (or rnapolis, or fr3d,rnapolis)."

    structure_ch = Channel.fromPath(input, checkIfExists: true)

    // T1 Conversion: any structure -> standardized mmCIF
    mmcif_ch = CONVERT(structure_ch)

    // T2 Annotate: one process per selected tool
    tools_ch = Channel.fromList(tools)
    raw_ch   = ANNOTATE( tools_ch.combine(mmcif_ch) )

    // T2 Parse: tool-native output -> standardized base-pairing mmCIF
    bp_ch = PARSE( raw_ch )

    bp_ch.view { tool, cif -> "[${tool}] annotated mmCIF -> ${cif}" }

    // T3 Validate: cross-tool comparison — per-tool stats always; the comparison is
    // skipped when only one tool is selected (gated in ws1-validate).
    VALIDATE( bp_ch.map { tool, cif -> cif }.collect() )

    // T3 Visualize (optional, --visualize): per-tool base-pair export — CSV/TSV tables
    // plus a forgi BulgeGraph (.bg) and a matplotlib render (.jpg).
    if( params.visualize )
        VISUALIZE( bp_ch )
}
