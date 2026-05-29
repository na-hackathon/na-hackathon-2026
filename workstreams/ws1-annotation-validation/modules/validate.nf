// T3 Validate — reparse the per-tool base-pair mmCIF (from PARSE) and emit
//   <tool>.csv  pair table (CSV)
//   <tool>.tsv  pair table (TSV)
//   <tool>.bg   forgi BulgeGraph built from cWW pairs
//   <tool>.jpg  matplotlib render of that BulgeGraph

process VALIDATE {
    tag   "${tool}"
    label 'validation'
    publishDir "${params.outdir}/validation", mode: 'copy'
    conda "${projectDir}/envs/parse.yml"

    input:
    tuple val(tool), path(bp_cif)

    output:
    tuple val(tool), path("${tool}.csv"), path("${tool}.tsv"), path("${tool}.bg"), path("${tool}.jpg")

    script:
    """
    python "${projectDir}/bin/export_basepairs.py" "${bp_cif}" "${tool}"
    """
}
