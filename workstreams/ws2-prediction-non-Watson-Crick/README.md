# Workstream 2 — RNA prediction & non-Watson-Crick base pairs

## Goals

- Interoperable representations for Watson-Crick + non-Watson-Crick base pairs
- Connect 2D representation to mmCIF
- Support benchmarking and case studies

## Where to start

- Open Issues labeled `workstream-2`
- Draft/iterate spec notes under this folder and `docs/`
## Test in iCn3D with R2DT, forna, and RNAcanvas

- Install [static-server](https://www.npmjs.com/package/static-server) and run "static-server -i module.html -o" to test your code. "import" and "export" do not work in "file://" protocol.
- Specify the ID of the structure in the URL, e.g., "localhost:9080/module.html?mmdbid=9CFN".
- Select the menu "Analysis > 2D Diagram > for Nucleotides", select the chain "9CFN_A", then click one of the button. For R2DT, you might want to try the chain 9 of PDB 1FFK. The eventual RNAcanvas might look like the following:
  <img width="398" height="401" alt="image" src="https://github.com/user-attachments/assets/f4c9a677-2b73-4b57-9816-e497c15c438b" />

