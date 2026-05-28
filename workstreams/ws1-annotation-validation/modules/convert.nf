// T1 Conversion — any input structure (.pdb/.cif) -> standardized mmCIF.
// Conversion team owns bin/ws1-convert and envs/convert.yml.

process CONVERT {
    tag   "${structure.name}"
    label 'conversion'
    publishDir "${params.outdir}/mmcif", mode: 'copy'
    conda "${projectDir}/envs/convert.yml"

    input:
    path structure

    output:
    path "${structure.baseName}.std.cif"

    script:
    """
    ws1-convert --input ${structure} --out ${structure.baseName}.std.cif
    """
}
