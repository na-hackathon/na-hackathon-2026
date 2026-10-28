// T2 Core (annotator) — run one base-pair annotator on the mmCIF.
// Calls bin/ws1-annotate, which dispatches per tool and runs the annotator from
// its pinned submodule under annotations_tools/ (passed via --src).
// Emits (tool, raw output, mmCIF, <tool>_provenance.json): the provenance file records the
// submodule commit, the command line and the sha256 of the mmCIF.

process ANNOTATE {
    tag   "${tool}:${mmcif.name}"
    label 'annotation'
    publishDir { "${params.outdir}/annotations/${tool}" }, mode: 'copy'
    conda "${projectDir}/envs/ws1.yml"

    input:
    tuple val(tool), path(mmcif)

    output:
    tuple val(tool), path("${tool}_raw.*"), path(mmcif), path("${tool}_provenance.json")

    script:
    """
    ws1-annotate --tool ${tool} --structure ${mmcif} --out ${tool}_raw --src ${projectDir}/annotations_tools \\
        --provenance ${tool}_provenance.json
    """
}
