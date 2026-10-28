// T3 Exploration
//   VALIDATE  : per-tool pair counts, and the cross-tool comparison of the base-pairing
//               mmCIFs when more than one tool is selected (gated in bin/ws1-validate).
//   VISUALIZE : optional per-tool export (params.visualize): CSV/TSV pair tables, and per
//               nucleic-acid chain a forgi BulgeGraph (.bg) and a matplotlib render (.jpg)
//               of the cWW secondary structure.

process VALIDATE {
    label 'validation'
    publishDir "${params.outdir}/validation", mode: 'copy'
    conda "${projectDir}/envs/ws1.yml"

    input:
    path bp_cifs

    output:
    path "validation_report.json"
    path "validation_report.tsv"

    script:
    """
    ws1-validate --annotations ${bp_cifs} --out validation_report
    """
}

process VISUALIZE {
    tag   "${tool}"
    label 'validation'
    publishDir "${params.outdir}/validation", mode: 'copy'
    conda "${projectDir}/envs/ws1.yml"

    input:
    tuple val(tool), path(bp_cif)

    output:
    tuple val(tool), path("${tool}.csv"), path("${tool}.tsv"), path("${tool}.bg"), path("${tool}.jpg")

    script:
    """
    python "${projectDir}/bin/export_basepairs.py" "${bp_cif}" "${tool}"
    """
}
