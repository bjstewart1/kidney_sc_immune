# kidney_sc_immune
Repository for scRNAseq study of human kidneys
Not a supported package.

## Anaylsis pipeline

This repository contains the following files which illustrate the methodology in our paper.
The functions and approaches described are assemmbled into a pipeline for the analysis of the data.

The following files pertain to this:
* Functions.R (functions used - mostly wrappers for convenience)
* SC_analysis_pipeline.R (single cell analysis pipeline)
* CellphoneDB_plotting.R (plotting pipeline for CellphoneDB)
* Pseudodepth_scores.R (approach to generate pseudodepth scores)

## Genesets
This folder contains genesets used from external sources for use in the analysis. Details in the README.md

## Data
This folder contains:
* Annotated final sce objects for fetal and mature kidney
* Microarray dataset (kidney depth study - Lindgren et al.)
* TCGA data (processed)
