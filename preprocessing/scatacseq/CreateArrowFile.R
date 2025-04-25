#!/usr/bin/ Rscript


####################
# Import libraries #
####################

library(argparse)
library(ArchR)

suppressPackageStartupMessages(library(ArchR))

##################################
# Define command line parameters #
##################################

parser <- ArgumentParser()

parser$add_argument("--data_dir", type = "character",
                    help = "Directory with scATAC-seq fragment files.")

parser$add_argument("--sample", type = "character",
                    help = "Sample name.")

parser$add_argument("--outdir", type = "character", 
		    help = "Output directory")

args <- parser$parse_args()
data_dir <- args$data_dir
sample <- args$sample 
outdir <- args$outdir

print("Selected arguments:\n")
print(paste0("Data dir: ", data_dir))
print(paste0("Sample name: ", sample))
print(paste0("Outdir: ", outdir))

################################
# Create ArrowFiles per sample #
################################

setwd(paste0(outdir, sample, "/"))

addArchRGenome("hg38")
addArchRThreads(threads = 8) 
library(parallel)

rhdf5::h5disableFileLocking()

ArrowFile <- createArrowFiles(
  inputFiles = paste0(data_dir, sample, '/fragments.tsv.gz'),
  sampleNames = sample,
  outputNames = sample,
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE, 
  force = TRUE
)

print("Created ArrowFile") 

##########################
# Compute doublet scores #
##########################

doubScores <- addDoubletScores(
    input = ArrowFile,
    k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
    knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
    LSIMethod = 1,
    outDir = paste0(outdir, sample, "/QC"),
    logFile = createLogFile("addDoubletScores", logDir = paste0(outdir, sample, '/QC'))
)




