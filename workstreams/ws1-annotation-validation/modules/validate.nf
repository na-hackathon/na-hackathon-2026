// T3 Exploration —
//   VALIDATE  : compare the per-tool base-pairing mmCIFs across tools. Per-tool pair
//               counts always; the cross-tool agreement comparison runs only when more
//               than one tool is selected (gated in bin/ws1-validate).
//   VISUALIZE : optional per-tool export (enabled by params.visualize) — CSV/TSV pair
//               tables plus a forgi BulgeGraph (.bg) and a matplotlib render (.jpg) of
//               the cWW secondary structure. From PR #86.

process VALIDATE {
    label 'validation'
    publishDir "${params.outdir}/validation", mode: 'copy'
    conda "${projectDir}/envs/parse.yml"

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
    conda "${projectDir}/envs/visualize.yml"

    input:
    tuple val(tool), path(bp_cif)

    output:
    tuple val(tool), path("${tool}.csv"), path("${tool}.tsv"), path("${tool}.bg"), path("${tool}.jpg")

    script:
    """
    python "${projectDir}/bin/export_basepairs.py" "${bp_cif}" "${tool}"
    """
}
