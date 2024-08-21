#!/usr/bin/env nextflow
  
params.cellbender = "/lustre/scratch127/cellgen/cellgeni/tickets/tic-1313/data19/starsolo/"
params.figures = "/nfs/team292/vl6/FetalReproductiveTract/RNA_QC_CellBender/figures/"
params.outdir = "/nfs/team292/vl6/FetalReproductiveTract/RNA_QC_CellBender/data/"

samples = Channel.from("HD_F_GON14896472")

project_dir = projectDir

log.info """\
         S C R N A S E Q - N F   P I P E L I N E
         ===================================
         cellbender directory: ${params.cellbender}
         figures directory: ${params.figures}
	 outdir: ${params.outdir}
         """
         .stripIndent()

process QC {

    input:
    val sample from samples 
    
    output:
    file ('*.h5ad') into adata_channel

    script:
    """
    python $project_dir/QC.py --sample $sample --cellbender ${params.cellbender} --figures ${params.figures}
    """
}

process PREPROCESSING {

    input:
    file adata from adata_channel 

    script:
    """
    python $project_dir/processing.py --anndata $adata --outdir ${params.outdir} 
    """
}

