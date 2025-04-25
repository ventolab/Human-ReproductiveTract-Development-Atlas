This directory contains the following scripts

## Per-library analyses 

- `SingleSample_ATAC_Preprocessing_Nextflow.nf`: Nextflow script that creates [ArchR](https://doi.org/10.1038/s41588-021-00790-6) Arrow files per sample
- `CreateArrowFile.R`: helper R script called from the Nextflow script
- `SingleSample_ATAC_ArchR.ipnyb`: R notebook to annotate cell types in each individual scATAC-seq library based on the matched view of the scRNA-seq dataset (e.g. if the scATAC-seq library is from an 8 PCW sample, we use the "early" <=10 PCW integrated scRNA-seq manifold to transfer cell type annotations) using [ArchR](https://doi.org/10.1038/s41588-021-00790-6)

## Integrations 

- `Integration_ATAC_Early_ArchR.ipnyb`: R notebook to integrate early (<=10 PCW) male and female samples using [ArchR](https://doi.org/10.1038/s41588-021-00790-6) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
- `Integration_ATAC_FemalesLate_ArchR.ipnyb`: R notebook to integrate late (>10 PCW) female samples using [ArchR](https://doi.org/10.1038/s41588-021-00790-6) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
- `Integration_ATAC_MalesLate_ArchR.ipnyb`: R notebook to integrate late (>10 PCW) male samples using [ArchR](https://doi.org/10.1038/s41588-021-00790-6) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
