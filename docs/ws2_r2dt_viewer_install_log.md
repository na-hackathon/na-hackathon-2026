
# Installation notes (Ubuntu)

## Clone R2DT repo

git clone https://github.com/r2dt-bio/R2DT.git

## Go to the R2DT folder

cd R2DT/

## Create output folder

mkdir output

## Install docker (if not yet)

sudo apt install docker.io

## Resolve "permission denied" issue for docker

sudo usermod -aG docker $USER

newgrp docker

## Run the 2d_3d viewer of R2DT (PDB id input)

docker run --rm -v $(pwd)/output:/rna/r2dt/output -w /rna/r2dt rnacentral/r2dt:pr-219 ./r2dt.py pdb_2d_3d 9RJA output/9RJA_2d3d

## Run the backend for the web results

python3 -m http.server -d ./output/9RJA_2d3d/viewer 8000

## Open the results in the browser

open http://localhost:8000/

## Copy the local mmCIF file to the output directory

cp [PATHTOFILE]/9CFN_clean.cif  ./output/9CFN_clean.cif

## Run the 2d_3d viewer of R2DT (local mmCIF file input)

docker run --rm -v $(pwd)/output:/rna/r2dt/output -w /rna/r2dt rnacentral/r2dt:pr-219 ./r2dt.py pdb_2d_3d /rna/r2dt/output/9CFN_clean.cif output/9CFN_2d3d
