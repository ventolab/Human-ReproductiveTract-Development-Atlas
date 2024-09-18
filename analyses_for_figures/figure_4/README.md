## Figure 4: Regulators of Müllerian and Wolffian duct mesenchymal patterning

This directory contains the following scripts 

- `1a_MüllerianAxis_Visium.ipynb`: python notebook that derives the computational representation of the differentiating human Müllerian duct based on landmarks annotated in the histology images of 10x Visium slides. 

- `1b_WolffianAxis_Visium.ipynb`: python notebook that derives the computational representation of the (partial) differentiating human Wolffian duct based on landmarks annotated in the histology images of 10x Visium slides. 
  
- `2a_MüllerianAxis_ISS.ipynb`: python notebook that derives the computational representation of the differentiating human Müllerian duct based on landmarks annotated in the virtual RGB images of In Situ Sequencing slides. It also performes the imputation of the position along the differentiating Müllerian duct in dissociated scRNA-seq data based on In Situ Sequencing data.  
  
- `3a_MüllerianAxis_Mesenchymal_Spatially_Variable_Genes_and_Spatially_Variable_Interactions.ipynb`: python notebook that intersects the mesenchymal spatially variable genes along the Müllerian rostro-caudal axis identified in 10x Visium and In Situ Sequencing data and computes spatially variable mesenchymal-epithelial interactions.

- `3b_WolffianAxis_Epithelial_Spatially_Variable_Genes.ipynb`: python notebook that prioritises the epithelial spatially variable genes along the partial rostro-caudal Wolffian duct axis identified in 10x Visium.

- `4a_MüllerianAxis_Epithelial_Spatially_Variable_Genes.ipynb`: python notebook that intersects the epithelial spatially variable genes along the Müllerian rostro-caudal axis identified in 10x Visium and In Situ Sequencing data.

- `4b_WolffianAxis_Mesenchymal_Spatially_Variable_Genes_and_Spatially_Variable_Interactions.ipynb`: python notebook that prioritises the mesenchymal spatially variable genes along the partial rostro-caudal Wolffian duct axis identified in 10x Visium and computes spatially variable mesenchymal-epithelial interactions.
