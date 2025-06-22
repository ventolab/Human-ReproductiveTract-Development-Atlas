This folder contains representative notebooks for stitching H&E images and corresponding anndata objects of consecutive 10x Visium slides:

- `ImageStitching_1_TransformationMatrices.ipynb`: apply transformation matrices to H&E images and to anndata objects using the JSON file with the spot coordinates to pixel correspondence for two consecutive sections of a representative second trimester female sample (17 PCW uterovaginal canal). Also save the overlapping spatial barcodes between the two sections being stitched.
- `ImageStitching_2_JoinSlides.ipynb`: create the stitched anndata object for two consecutive sections of a representative second trimester female sample (17 PCW uterovaginal canal)
- `ImageStitching_3_Clustering.ipynb`: compute scale factors between overlapping spatial barcodes and normalise gene expression for the stitched anndata object based on the scale factors. Cluster the resulting gene expression matrix and plot genes for the stitched anndata object for two consecutive sections of a representative second trimester female sample (17 PCW uterovaginal canal)

