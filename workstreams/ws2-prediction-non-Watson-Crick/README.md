# Workstream 2 — RNA prediction & non-Watson-Crick base pairs

## Goals

- Interoperable representations for Watson-Crick + non-Watson-Crick base pairs
- Connect 2D representation to mmCIF
- Support benchmarking and case studies

## Where to start

- Open Issues labeled `workstream-2`
- Draft/iterate spec notes under this folder and `docs/`
## Test in iCn3D with R2DT, forna, and RNAcanvas

- Install [static-server](https://www.npmjs.com/package/static-server) and run "static-server -i rna2d.html -o" to test your code. "import" and "export" do not work in "file://" protocol.
- Specify the ID of the structure in the URL, e.g., "http://localhost:9080/rna2d.html?mmdbid=9CFN". You can also load mmcif file in the same directory with this URL: "http://localhost:9080/rna2d.html?urltype=mmcif&urlname=http://localhost:9080/9hrf.cif".
- Select the menu "Analysis > 2D Diagram > for Nucleotides", select the chain "9CFN_A", then click one of the button. For R2DT, you might want to try the chain 9 of PDB 1FFK.
## 2D diagram comparison for chain 0 of PDB 1FFK
1. R2DT: <img width="661" height="762" alt="image" src="https://github.com/user-attachments/assets/35077bf6-3c7d-469a-af99-c9cbd871153e" />
2. forna: <img width="882" height="788" alt="image" src="https://github.com/user-attachments/assets/2071a002-c1b7-437b-a574-215467849e96" />
3. RNAcanvas: <img width="632" height="690" alt="image" src="https://github.com/user-attachments/assets/7a37bb9e-37ec-4188-a020-727e6be209ca" />
Plan: R2DT is the best. RNAcanvas is not bad. forna is too messy.





