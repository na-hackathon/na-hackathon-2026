#!/usr/bin/env nextflow

// At the top of your main.nf file
log.info "The project directory is: ${projectDir}"

process RNApolisAnnotate {
    tag "${structure.name}"
    label 'annotation'
    publishDir "${params.outdir}/annotation/rnapolis", mode: 'copy'
    conda "env.yml"
    
    input:
    path structure
    
    output:
    path "${structure.simpleName}.csv", emit: base_pairs_csv
    path "${structure.simpleName}.json", emit: base_pairs_json

    script:
    // Using a Nextflow variable to keep the bash block clean
    def prefix = structure.simpleName
    """
    annotator --csv "${prefix}.csv" --json "${prefix}.json" --extended "${structure}"
    """
}

workflow {    
    // Ensure params.input is defined or fallback to avoid a null error
    def input_path = params.input ?: 'data/*.cif.gz'
    
    structure_ch = Channel.fromPath(input_path, checkIfExists: true)
    raw_ch       = RNApolisAnnotate( structure_ch )
}
