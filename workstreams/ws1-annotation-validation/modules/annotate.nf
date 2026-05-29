// T2 Core (annotator) — run one base-pair annotator on the mmCIF.
// Calls bin/ws1-annotate, which dispatches per tool and runs the annotator from
// its pinned submodule under annotations_tools/ (passed via --src).

process ANNOTATE {
    tag   "${tool}:${mmcif.name}"
    label 'annotation'
    publishDir { "${params.outdir}/annotations/${tool}" }, mode: 'copy'
    conda { "${projectDir}/envs/${tool}.yml" }

    input:
    tuple val(tool), path(mmcif)

    output:
    tuple val(tool), path("${tool}_raw.*"), path(mmcif)

    script:
    """
    ws1-annotate --tool ${tool} --structure ${mmcif} --out ${tool}_raw --src ${projectDir}/annotations_tools
    """
}
