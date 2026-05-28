// T3 Exploration — compare the standardized base-pairing mmCIFs across tools and
// emit a validation report. Validation team owns bin/ws1-validate.

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
