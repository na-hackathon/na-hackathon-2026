


## We

Maciej Antczak, Eugene Baulin, Zubaida Khatoon, Sophie Maclure, Patrik Marek, Amrit Patra, Anton Petrov, Jim Procter, Agnieszka Rybarczyk, Joanna Sarzynska, David Sehnal, Kunal Shewani, Marta Szachniuk, Jiyao Wang

## We started with an overview of existing solutions

https://docs.google.com/document/d/1z3ibb5wq2gQKVdQtDDs0Qbp1jKlpfLVuBvKr1a_D20k/edit?tab=t.0


## We ended up with a specification of a 2D-3D interactive visualization

https://github.com/febos/na-hackathon-2026/blob/main/docs/ws2_visualization_spec.md

## Direction 1 - icn3D + forna/RNAcanvas (lead: Jiyao Wang)

### GitHub location

https://github.com/na-hackathon/na-hackathon-2026/tree/main/workstreams/ws2-prediction-non-Watson-Crick

### How to run locally

static-server -i rna2d.html -o

http://localhost:9080/rna2d.html?mmdbid=9CFN

http://localhost:9080/rna2d.html?urltype=mmcif&urlname=http://localhost:9080/9hrf.cif

## Direction 2 - layered dot-bracket (lead: Amrit Patra)

### GitHub location

https://github.com/na-hackathon/na-hackathon-2026/tree/main/workstreams/ws2-prediction-non-Watson-Crick/layered-bp-notation

## Direction 3 - Mol* + R2DT + dot-bracket (lead: Anton Petrov)

### Web demo

https://2d-3d-viewer.na-hackathon-2026.pages.dev/workstream1/

### GitHub location

https://github.com/r2dt-bio/R2DT/blob/2d-3d-viewer/docs/pdb-2d-3d-viewer.md#running-via-docker

### How to run locally with regular mmCIF

docker run --rm -v $(pwd)/output:/rna/r2dt/output -v $(pwd):/rna/r2dt -w /rna/r2dt rnacentral/r2dt:pr-219 ./r2dt.py pdb_2d_3d /rna/r2dt/output/9CFN_clean.cif  output/9CFN_2d3d

python3 -m http.server -d ./output/9CFN_2d3d/viewer 8000

http://0.0.0.0:8000/

### How to run locally with extended mmCIF

docker run --rm -v $(pwd)/output:/rna/r2dt/output -v $(pwd):/rna/r2dt -w /rna/r2dt rnacentral/r2dt:pr-219 ./r2dt.py pdb_2d_3d --basepairs cif /rna/r2dt/output/9CFN_dnatco.cif  output/9CFN_2d3d

python3 -m http.server -d ./output/9CFN_2d3d/viewer 8000

http://0.0.0.0:8000/

## Direction 4 - optimized zoom on base pairs in Mol* (lead: David Sehnal)

### GitHub location of demo video

https://github.com/na-hackathon/na-hackathon-2026/blob/main/docs/molstar_basepair_zoom.mov


