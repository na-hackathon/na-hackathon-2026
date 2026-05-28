process GEMMI_CONVERTPDB2MMCIF {
    tag   "${structure.name}"
    label 'conversion'
    publishDir "${params.outdir}/mmcif", mode: 'copy'
    conda "${projectDir}/envs/convert.yml"

    input:
    path structure

    output:
    path "${structure.baseName}.std.cif"

    script:
    if (structure.extension == 'pdb')
        """
        gemmi convert --from=pdb --to=mmcif ${structure} "${structure.baseName}.std.cif"
        """
    else
        """
        cat ${structure} > ${structure.baseName}.std.cif
        """

}