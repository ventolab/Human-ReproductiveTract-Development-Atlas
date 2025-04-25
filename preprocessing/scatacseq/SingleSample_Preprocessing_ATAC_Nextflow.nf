#!/usr/bin/env nextflow

params.data_dir = "/nfs/team292/vl6/FetalReproductiveTract/ATAC_QC/data/"
params.outdir = "/nfs/team292/vl6/FetalReproductiveTract/ATAC_QC/ArchR/"

samples = Channel.from( "HD_F_GON15261138"
)

project_dir = projectDir

log.info """\
         S C A T A C S E Q - N F   P I P E L I N E
         =========================================
         data directory: ${params.data_dir}
         outdir: ${params.outdir}
         """
         .stripIndent()


process ARROWFILE {

    input:
    val sample from samples

    script:
    """
    Rscript $project_dir/CreateArrowFile.R --sample $sample --data_dir ${params.data_dir} --outdir ${params.outdir}
    """
}
