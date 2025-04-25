This directory contains the following scripts

## Per-library analysis 

- `scatacseq_preprocessing_nextflow.nf`
- `SingleSample_ATAC_ArchR.ipnyb`: R notebook to annotate cell types in each individual scATAC-seq library based on the matched view of the scRNA-seq dataset (e.g. if the scATAC-seq library is from an 8 PCW sample, we use the "early" <=10 PCW integrated scRNA-seq manifold to transfer cell type annotations)
