This directory contains the following scripts

## Per-library analyses and integration of scATAC-seq data

- `scatacseq_preprocessing_nextflow.nf`
- `SingleSample_ATAC_ArchR.ipnyb`: R notebook to annotate cell types in each individual scATAC-seq library based on the matched view of the scRNA-seq dataset (e.g. if the scATAC-seq library is from an 8 PCW sample, we use the "early" <=10 PCW integrated scRNA-seq manifold to transfer cell type annotations) using [ArchR]([10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6))
- `Integration_ATAC_Early_ArchR.ipnyb`: R notebook to integrate early (<=10 PCW) male and female samples using [ArchR]([10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6)) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
- `Integration_ATAC_FemalesLate_ArchR.ipnyb`: R notebook to integrate late (>10 PCW) female samples using [ArchR]([10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6)) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
- `Integration_ATAC_MalesLate_ArchR.ipnyb`: R notebook to integrate late (>10 PCW) male samples using [ArchR]([10.1038/s41588-021-00790-6](https://doi.org/10.1038/s41588-021-00790-6)) and [Mutual Nearest Neighbours](https://doi.org/10.1038/nbt.4091)
