#!/usr/bin/env nextflow

// End-to-end annotation pipeline:
//   structure (.pdb|.cif) -> T1 CONVERT (MaxIT) -> T2 ANNOTATE (fr3d|rnapolis)
//   -> T2 PARSE -> standardized base-pairing mmCIF (per tool)
//   -> T3 VALIDATE (cross-tool comparison) + optional VISUALIZE (--visualize: per-tool graphs)
//   RUN_MANIFEST records versions, commits and checksums of the run.
//
// Usage:
//   nextflow run main.nf <structure.pdb|.cif>                              # both tools, default
//   nextflow run main.nf --input <file> --annotators fr3d                  # single tool
//   nextflow run main.nf --input <file> --annotators fr3d,rnapolis         # both tools + validate
//
// The per-tool annotated mmCIF is published to ${params.outdir}/basepairs/${tool}_basepairs.cif.

include { CONVERT      } from './modules/convert.nf'
include { ANNOTATE     } from './modules/annotate.nf'
include { PARSE        } from './modules/parse.nf'
include { VALIDATE; VISUALIZE } from './modules/validate.nf'
include { RUN_MANIFEST } from './modules/manifest.nf'

// "a=1, b=2" for a map of counts
def formatCounts(counts) {
    return counts.collect { name, n -> "${name}=${n}" }.join(', ')
}

// Console summary of validation_report.json; the full table is validation_report.tsv.
def validationSummary(report) {
    def r = new groovy.json.JsonSlurper().parse(report.toFile())
    def lines = [ "tools:                ${r.tools.join(', ')}",
                  "pairs per tool:       ${formatCounts(r.pair_counts)}" ]
    if( r.total_distinct_pairs != null ) {
        lines << "total distinct pairs: ${r.total_distinct_pairs}"
        lines << "agreed by all tools:  ${r.agreed_by_all_tools} (discordant ${r.discordant_family}; only in one tool: ${formatCounts(r.only_in)})"
    }
    lines << "full report:          ${params.outdir}/validation/validation_report.tsv"
    return lines.join('\n')
}

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
    PARSE( raw_ch )
    bp_ch = PARSE.out.basepairs

    bp_ch.view { tool, cif -> "[${tool}] annotated mmCIF -> ${cif}" }

    // T3 Validate: cross-tool comparison — per-tool stats always; the comparison is
    // skipped when only one tool is selected (gated in ws1-validate).
    VALIDATE( bp_ch.map { tool, cif -> cif }.collect() )
    VALIDATE.out.json.view { report -> validationSummary(report) }

    // T3 Visualize (optional, --visualize): per-tool base-pair export — CSV/TSV tables
    // plus, per nucleic-acid chain, a forgi BulgeGraph (.bg) and a matplotlib render (.jpg).
    if( params.visualize )
        VISUALIZE( bp_ch )

    // Run manifest: Nextflow facts of this run + the provenance of every annotator
    run_info = Channel.of( groovy.json.JsonOutput.toJson([
        nextflow    : workflow.nextflow.version.toString(),
        command_line: workflow.commandLine,
        run_name    : workflow.runName,
        profile     : workflow.profile,
        start       : workflow.start.toString(),
        params      : params.collectEntries { k, v -> [k, v?.toString()] }
    ]) ).collectFile(name: 'run_info.json')
    RUN_MANIFEST( structure_ch, raw_ch.map { tool, raw, mmcif, provenance -> provenance }.collect(), run_info )
}
