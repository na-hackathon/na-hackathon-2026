// Run manifest: pipeline commit, Nextflow facts, input checksum, annotator provenance and the
// package versions of the process env, written to <outdir>/run_manifest.json (bin/ws1-manifest).

process RUN_MANIFEST {
    label 'manifest'
    publishDir "${params.outdir}", mode: 'copy'
    conda "${projectDir}/envs/ws1.yml"

    input:
    path structure          // the input of the run
    path provenance         // <tool>_provenance.json of every annotator
    path run_info           // run_info.json, written by main.nf

    output:
    path "run_manifest.json"

    script:
    """
    ws1-manifest --input ${structure} --run-info ${run_info} --project-dir ${projectDir} \\
        --provenance ${provenance} --out run_manifest.json
    """
}
