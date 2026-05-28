

# Specification for a prototype of interactive 2D/3D visualizations

- based on the discussion of existing solutions, 27.05.26 (https://docs.google.com/document/d/1z3ibb5wq2gQKVdQtDDs0Qbp1jKlpfLVuBvKr1a_D20k/edit)

## Required (minimum viable version)

1) 2D viewer + 3D viewer + interactivity between the two (minimum intractivity = clickability on residue/base pair; optional extensions = e.g., show all base pairs of the residue on click)

2) Reproducibility (it can be installed and run locally)

3) An option to toggle canonical/non-canonical base pairs (canonical-only or all)

4) Show residue identifier on mouse hover in both panels (as pop-up)

## Optional (good to have)

1) 1D viewer (Rchie/dot-bracket) + 2D viewer + 3D viewer + interactivity across them

2) Feature base pair symbols (square, circle, triangle) in the 2D viewer

3) Feature a base pair list as an additional panel (possibly also interactive)

4) Toggle nested / non-nested base pairs

5) Select a subset of the input mmCIF (one chain, one motif, etc.) - before the visualization (apply to all the panels) or after (apply to some of the panels)?

6) Handle "near" base pairs and other "unusual" cases. For example, handle an extended FR3D annotation (see 9CFN_bgsu_fr3d_all.tsv)


## Current directions

 - R2DT + Mol* (by Anton Petrov)

 - icn3D2 + forna (by Jiyao Wang)
 
 - Jalview (by Jim Procter)
 
 Additional example structures (large ones? modified residues? DNA?) and feedback are very welcome!
