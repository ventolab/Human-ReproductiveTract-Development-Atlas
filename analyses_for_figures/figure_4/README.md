## Figure 4: Regulators of Müllerian and Wolffian duct mesenchymal patterning & Figure 5: Cell-cell communication along the Müllerian and Wolffian duct niches

This directory contains the following scripts 

**Müllerian and Wolffian rostro-caudal axis construction**

- `1a_MüllerianAxis_Visium.ipynb`: python notebook that derives the computational representation of the differentiating human Müllerian duct based on landmarks annotated in the histology images of 10x Visium slides. 

- `1b_WolffianAxis_Visium.ipynb`: python notebook that derives the computational representation of the (partial) differentiating human Wolffian duct based on landmarks annotated in the histology images of 10x Visium slides. 
  
- `2a_MüllerianAxis_ISS.ipynb`: python notebook that derives the computational representation of the differentiating human Müllerian duct based on landmarks annotated in the virtual RGB images of In Situ Sequencing slides. It also performes the imputation of the position along the differentiating Müllerian duct in dissociated scRNA-seq data based on In Situ Sequencing data.

**Sub-analyses of Müllerian and Wolffian-derived mesenchymal and epithelial cells (scRNA-seq data)**
- `3a_Müllerian_Duct_Differentiation_Mesenchyme.ipynb`:
  
- `3b_Wolffian_Duct_Differentiation_Mesenchyme.ipynb`:
  
- `4a_Müllerian_Duct_Differentiation_Epithelium.ipynb`:
  
- `4b_Wolffian_Duct_Differentiation_Epithelium.ipynb`:
  
- `5a_Müllerian_Duct_Differentiation_Mesenchyme+Epithelium.ipynb`:
  
- `5b_Wolffian_Duct_Differentiation_Mesenchyme+Epithelium.ipynb`: 

**Continuous modelling of gene expression along the Müllerian and Wolffian rostro-caudal axes**
  
- `6a_MüllerianAxis_Mesenchymal_Spots_TradeSeq.ipynb`:  R notebook that models spatially-variable genes in mesenchymal spots (from 10x Visium) along the measured Müllerian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)

- `6b_WolffianAxis_Mesenchymal_Spots_TradeSeq.ipynb`: R notebook that models spatially-variable genes in mesenchymal spots (from 10x Visium) along the measured Wolffian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)

- `7a_MüllerianAxis_Mesenchymal_Cells_TradeSeq.ipynb`: R notebook that models spatially-variable genes in mesenchymal cells (from scRNA-seq) along the imputed Müllerian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)
  
- `8a_MüllerianAxis_Epithelial_Spots_TradeSeq.ipynb`: R notebook that models spatially-variable genes in epithelial spots (from 10x Visium) along the measured Müllerian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)

- `8b_WolffianAxis_Epithelial_Spots_TradeSeq.ipynb`: R notebook that models spatially-variable genes in epithelial spots (from 10x Visium) along the measured Wolffian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)

- `9a_MüllerianAxis_Epithelial_Cells_TradeSeq.ipynb`: R notebook that models spatially-variable genes in epithelial cells (from scRNA-seq) along the imputed Müllerian rostro-caudal axis using [TradeSeq](https://www.nature.com/articles/s41467-020-14766-3)

**Prioritisation of spatially-variable TFs and cell-cell communicational events between mesenchyme and epithelium**

- `10a_MüllerianAxis_Mesenchymal_Spatially_Variable_Genes_and_Spatially_Variable_Interactions.ipynb`: python notebook that intersects the mesenchymal spatially variable genes along the Müllerian rostro-caudal axis identified in 10x Visium and In Situ Sequencing data and computes spatially variable mesenchymal-epithelial interactions.

- `10b_WolffianAxis_Mesenchymal_Spatially_Variable_Genes_and_Spatially_Variable_Interactions.ipynb`: python notebook that prioritises the mesenchymal spatially variable genes along the partial rostro-caudal Wolffian duct axis identified in 10x Visium and computes spatially variable mesenchymal-epithelial interactions.

- `11a_MüllerianAxis_Epithelial_Spatially_Variable_Genes.ipynb`: python notebook that intersects the epithelial spatially variable genes along the Müllerian rostro-caudal axis identified in 10x Visium and In Situ Sequencing data.

- `11b_WolffianAxis_Epithelial_Spatially_Variable_Genes.ipynb`: python notebook that prioritises the epithelial spatially variable genes along the partial rostro-caudal Wolffian duct axis identified in 10x Visium.
