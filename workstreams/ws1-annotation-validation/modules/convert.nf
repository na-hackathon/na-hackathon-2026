// T1 Conversion — any input structure (.pdb/.cif) -> standardized mmCIF.
// Conversion team owns bin/ws1-convert and envs/convert.yml.

process CONVERT {
    tag   "${structure.name}"
    label 'conversion'
    publishDir "${params.outdir}/mmcif", mode: 'copy'
    //conda "${projectDir}/envs/convert.yml"
    container "tzok/maxit"

    input:
    path structure

    output:
    path "${structure.baseName}.std.cif"

    script:
    """
    maxit -input ${structure} -output ${structure.baseName}.std.cif -o 1
    """
}


// python ${projectDir}/bin/pdb_to_mmcif.py -o ${structure.baseName}.std.cif ${structure}