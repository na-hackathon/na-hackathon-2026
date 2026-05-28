// T1 Conversion — any input structure (.pdb/.cif) -> standardized mmCIF.
// Conversion team owns bin/ws1-convert and envs/convert.yml.

process CONVERT {
    tag   "${structure.name}"
    label 'conversion'
    publishDir "${params.outdir}/mmcif", mode: 'copy'
    conda "${projectDir}/envs/convert.yml"
    //this defines a docker container in which the maxit can run on it's own, but that stops the python wrapping from working
    //container "tzok/maxit"

    input:
    path structure

    output:
    path "${structure.baseName}.std.cif"

    script:
    """
    python ${projectDir}/bin/convert_to_mmcif_with_maxit.py -o ${structure.baseName}.std.cif ${structure}
    """
    //if you want to use this process but not the python wrapped version, enable docker and swap out the script above with what's below
    //if (structure.extension == 'pdb')
    //    """
    //    maxit -input ${structure} -output ${structure.baseName}.std.cif -o 1
    //    """
    //else
    //    """
    //    maxit -input ${structure} -output ${structure.baseName}.std.cif -o 8
    //    """
}
