// T2 Core (annotator) — run one base-pair annotator on the mmCIF.
// Generic over the tool name: each tool supplies bin/ws1-annotate behaviour and
// its own envs/<tool>.yml. Add a tool by dropping in those two files and listing
// the name in --annotators. The mmCIF is carried through for the parser stage.

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
    ws1-annotate --tool ${tool} --structure ${mmcif} --out ${tool}_raw
    """
}
