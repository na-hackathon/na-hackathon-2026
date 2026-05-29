// T2 Core (parser) — tool-native annotation output -> standardized base-pairing
// mmCIF (NDB extended categories). Calls bin/ws1-parse, which dispatches per tool.

process PARSE {
    tag   "${tool}"
    label 'parsing'
    publishDir "${params.outdir}/basepairs", mode: 'copy'
    conda "${projectDir}/envs/parse.yml"

    input:
    tuple val(tool), path(raw), path(mmcif)

    output:
    tuple val(tool), path("${tool}_basepairs.cif")

    script:
    """
    ws1-parse --tool ${tool} --raw ${raw} --structure ${mmcif} --out ${tool}_basepairs.cif
    """
}
