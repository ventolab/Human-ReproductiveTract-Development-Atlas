This repository hosts the code and data processing pipelines used for the analysis of single-cell RNA sequencing (scRNA-seq), single-cell ATAC sequencing (scATAC-seq), and spatially-resolved transcriptomics by means of 10x Visium and In Situ Sequencing (ISS) data in the study of the developing human reproductive tract. The repository is structured to facilitate the reproduction of the results presented in our paper, including all figures and key analyses.

## Repository Structure
The repository is organized into the following main directories:

1. `preprocessing/`
This folder contains the code for preprocessing and initial analysis of the raw sequencing data. The preprocessing is divided into the following subdirectories:

- `scrnaseq/`: Scripts and workflows for processing single-cell RNA sequencing data.
- `scatacseq/`: Scripts and workflows for processing single-cell ATAC sequencing data.
- `10xVisium/`: Code for processing spatial transcriptomics data obtained using the 10x Visium platform. This includes the annotation of anatomical and histological structures in the H&E images.
- `ISS/`: Workflows and scripts for processing In Situ Sequencing data. This includes the annotation of anatomical structures in the virtual H&E images.

2. `analyses_for_figures/` (`figure_1/` to `figure_7/`)
These folders contain the analysis code required to reproduce each figure and the associated results presented in the manuscript. Below is a brief description of each figure directory, where you can add more detailed information about the specific analyses:

- `figure_1_atlas/`: Spatiotemporal atlas of human reproductive tract development 

- `figure_2_Müllerian_emergence/`: Ontology, migration and regression of the Müllerian ducts

- `figures_3_4_Müllerian_Wolffian_regionalisation/`: Regulators of Müllerian and Wolffian duct mesenchymal patterning & cell-cell communication along the Müllerian and Wolffian duct niches

- `figure_5_intra_organ_regionalisation/`: Regionalisation of the fallopian tube and epididymis occurs during fetal development

- `figure_6_external_genitalia/`: Sexual dimorphism in the genital tubercle

- `figure_7_organoids_disruption/`: Mapping potential disruptions to reproductive tract development & in vitro validation of estrogen-mimicking endocrine disrupting chemicals 
  
Each figure directory is self-contained and includes:

- Data preparation: Code for loading and preparing the processed data for analysis.
- Analysis scripts: Scripts that perform the necessary computations and statistical analyses.
- Figure generation: Code that generates the figures as they appear in the manuscript, including any supplementary materials.

