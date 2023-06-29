# singlecell_humanmelanoma
Code for scRNAseq of three syngeneic human melanoma cell lines (two replicate per cell line) to study mechanisms of resistance to the BRAF inhibitor vemurafenib

Raw and processed data can be downloaded in the following link:
https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-12412

- 0.1.readData_rep1.r and 0.2.readData_rep2.r for loading replicate 1 and replicate 2 data into R.
- 1.mergeData.r for merging the two replicates into a Seurat object.
- 2.filter_mithocondrial.r for filtering low quality cells and ribosomal genes.
- 3.analysis_VR_VRRANO.r for downstream analysis of VR vs VR-RANO conditions.
- 4.analysis_parental_VR_VRRANO.r for downstream analysis of Parental vs VR vs VR-RANO conditions.

- MDA_1.1.14.tar.gz, functions.r and utils.R contain in-house functions that are used in the main code. 
