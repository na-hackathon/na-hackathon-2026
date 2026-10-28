// T2 Core (parser) — tool-native annotation output -> standardized base-pairing
// mmCIF (NDB extended categories). Calls bin/ws1-parse, which dispatches per tool.
// Besides <tool>_basepairs.cif it emits <tool>_basepairs.skipped.tsv: the input pairs the
// exporter did not write, with the reason.

process PARSE {
    tag   "${tool}"
    label 'parsing'
    publishDir "${params.outdir}/basepairs", mode: 'copy'
    conda "${projectDir}/envs/ws1.yml"

    input:
    tuple val(tool), path(raw), path(mmcif)

    output:
    tuple val(tool), path("${tool}_basepairs.cif"), emit: basepairs
    path "${tool}_basepairs.skipped.tsv",           emit: skipped

    script:
    """
    ws1-parse --tool ${tool} --raw ${raw} --structure ${mmcif} --out ${tool}_basepairs.cif
    """
}
