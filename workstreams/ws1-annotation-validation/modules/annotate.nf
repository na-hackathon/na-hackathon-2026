// T2 Core (annotator) — run one base-pair annotator on the mmCIF.
// Dispatches per tool, invoking the source trees bundled under annotations_tools/:
//   fr3d     : annotations_tools/fr3d-python/fr3d/classifiers/NA_pairwise_interactions.py
//              -> rename the produced *_basepair.txt to ${tool}_raw.txt
//   rnapolis : annotations_tools/rnapolis-py/src/rnapolis/annotator.py -j <out.json> -e
// The mmCIF is carried through for the parser stage.

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
    def fr3d_src     = "${projectDir}/annotations_tools/fr3d-python"
    def rnapolis_src = "${projectDir}/annotations_tools/rnapolis-py/src"
    if( tool == 'fr3d' )
        """
        export PYTHONPATH="${fr3d_src}:\${PYTHONPATH:-}"
        python "${fr3d_src}/fr3d/classifiers/NA_pairwise_interactions.py" -i . ${mmcif} -o ${tool}_raw
        mv \$(find . -name '*_basepair.txt' | head -n 1) ${tool}_raw.txt
        # Having fixed name like the _basepair.txt is not super good but i don't have a lot of better ideas to fix the things.
        """
    else if( tool == 'rnapolis' )
        """
        export PYTHONPATH="${rnapolis_src}:\${PYTHONPATH:-}"
        python "${rnapolis_src}/rnapolis/annotator.py" ${mmcif} -j ${tool}_raw.json -e
        """
    else
        error "ANNOTATE: no annotator wired for tool '${tool}'. Add a branch in modules/annotate.nf."
}
