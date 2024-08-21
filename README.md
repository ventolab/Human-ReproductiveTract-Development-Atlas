This repository hosts the code and data processing pipelines used for the analysis of single-cell RNA sequencing (scRNA-seq), single-cell ATAC sequencing (scATAC-seq), and spatially-resolved transcriptomics by means of 10x Visium and In Situ Sequencing (ISS) data in the study of the developing human reproductive tract. The repository is structured to facilitate the reproduction of the results presented in our paper, including all figures and key analyses.

## Repository Structure
The repository is organized into the following main directories:

1. `preprocessing/`
This folder contains the code for preprocessing and initial analysis of the raw sequencing data. The preprocessing is divided into the following subdirectories:

- `scRNA-seq/`: Scripts and workflows for processing single-cell RNA sequencing data.
- `scATAC-seq/`: Scripts and workflows for processing single-cell ATAC sequencing data.
- `10xVisium/`: Code for processing spatial transcriptomics data obtained using the 10x Visium platform. This includes the annotation of anatomical and histological structures in the H&E images.
- `ISS/`: Workflows and scripts for processing In Situ Sequencing data. This includes the annotation of anatomical structures in the virtual H&E images.

2. `analyses_for_figures/` (`figure_1/` to `figure_7/`)
These folders contain the analysis code required to reproduce each figure and the associated results presented in the manuscript. Below is a brief description of each figure directory, where you can add more detailed information about the specific analyses:

- `figure_1/`: Spatiotemporal atlas of human reproductive tract development 

- `figure_2/`: Sexual dimorphism in the genital tubercle

- `figure_3/`: Ontology, migration and regression of the Müllerian ducts

- `figure_4/`: Regulators of Müllerian and Wolffian duct mesenchymal patterning

- `figure_5/`: Cell-cell communication along the Müllerian and Wolffian duct niches

- `figure_6/`: Regionalisation of the fallopian tube and epididymis occurs during fetal development

- `figure_7/`: Mapping potential disruptions to reproductive tract development
  
Each figure directory is self-contained and includes:

- Data preparation: Code for loading and preparing the processed data for analysis.
- Analysis scripts: Scripts that perform the necessary computations and statistical analyses.
- Figure generation: Code that generates the figures as they appear in the manuscript, including any supplementary materials.

