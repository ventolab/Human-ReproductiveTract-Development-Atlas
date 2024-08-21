This directory contains the following scripts 

## Per-library analysis 
- `scrnaseq_preprocessing_nextflow.nf`: nextflow pipeline with two processes to enable parallelisation and reproducibility. The first process (`QC.py`) computes and plots quality control metrics per library, while the second process (`processing.py`) carries out the preprocessing steps until the construction of the 2-D UMAP embedding. Each processed is applied to a given scRNA-seq library independently and it takes as input the output of STARSolo + CellBender.
- `QC.py`: script called by process one in the nextflow pipeline
- `processing.py`: script called by process two in the nextflow pipeline
- `reptract_utils.py`: handy functions for computing several of the analyses on scRNA-seq data
- `reptract_genes.py`: dictionary of the major marker genes per cell types in the reproductive tract that are used to guide cell type annotations
- `scrnaseq_annotation_per_sample.ipynb`: jupyter notebook to perform the preliminary, per-library cell type annotation based on marker genes curated from the literature 
