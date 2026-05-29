// T1 Conversion — standardize an input structure (.pdb/.cif) to mmCIF.
// Calls bin/ws1-convert. params.converter selects the converter and its env:
//   maxit  (default)    — RCSB MAXIT via conda (bioconda::maxit); pdb -> -o 1, cif -> -o 8 (#76)
//   maxit-docker        — same maxit from the tzok/maxit image (fallback)
//   gemmi  (DEPRECATED) — conda fallback
// maxit/gemmi come from envs/convert.yml; maxit-docker runs in the tzok/maxit container.

process CONVERT {
    tag   "${structure.name}"
    label 'conversion'
    publishDir "${params.outdir}/mmcif", mode: 'copy'
    conda     { params.converter == 'maxit-docker' ? null : "${projectDir}/envs/convert.yml" }
    container { params.converter == 'maxit-docker' ? 'tzok/maxit' : null }

    input:
    path structure

    output:
    path "${structure.baseName}.std.cif"

    script:
    """
    ws1-convert --input ${structure} --out ${structure.baseName}.std.cif --converter ${params.converter}
    """
}
