// T2 Core (parser) — tool-native annotation output -> standardized base-pairing
// mmCIF (NDB extended categories). Dispatches per tool to the conversion
// scripts in read_write_mmcif/.
//   fr3d     : raw TSV -> JSON -> NDB base-pairing mmCIF
//   rnapolis : raw JSON -> NDB base-pairing mmCIF

process PARSE {
    tag   "${tool}"
    label 'parsing'
    publishDir "${params.outdir}/basepairs", mode: 'copy'
    conda "${projectDir}/envs/parse.yml"

    input:
    tuple val(tool), path(raw), path(mmcif)

    output:
    tuple val(tool), path("${tool}_basepairs.cif")

    script:
    def scripts = "${projectDir}/read_write_mmcif"
    if( tool == 'fr3d' )
        """
        export PYTHONPATH="${scripts}:\${PYTHONPATH:-}"
        python3 "${scripts}/fr3d_basepair_to_json.py" "${raw}" basepairs.json
        python3 "${scripts}/json_fr3d_mmcif.py" basepairs.json "${mmcif}" "${tool}_basepairs.cif"
        """
    else if( tool == 'rnapolis' )
        """
        export PYTHONPATH="${scripts}:\${PYTHONPATH:-}"
        python3 "${scripts}/json_rnapolis_mmcif.py" "${raw}" "${mmcif}" "${tool}_basepairs.cif"
        """
    else
        error "PARSE: no parser wired for tool '${tool}'. Add a branch in modules/parse.nf or a script in read_write_mmcif/."
}
