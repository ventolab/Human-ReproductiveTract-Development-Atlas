#!/usr/bin/env python3

####################
# Import libraries #
####################

import argparse
import scrublet as scr
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy
sc.settings.verbosity = 3
sc.logging.print_versions()
sys.executable
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
##################################
# Import reptract utils function #
##################################

utils_dir = '/nfs/team292/vl6/RepTract/'
sys.path.append(utils_dir)

import reptract_utils

##################################
# Define command line parameters #
##################################

parser = argparse.ArgumentParser()

# anndata
parser.add_argument('--anndata', nargs='+',  type=str, help='Path to anndata object created by process QC of the Nextflow pipeline.', required = True)

# scrublet
parser.add_argument('--scrublet', nargs='+',  type=str, help='Path to where you want to save scrublet scores.', required = False)

# starsolo 
parser.add_argument('--starsolo', type=str, nargs='+', help='Path to STARsolo output.',required=False)

# outdir
parser.add_argument('--outdir', type=str, nargs='+', help='Path to directory where you want to save the resulting h5ad object.',required=True)

args = parser.parse_args()

####################
# Print parameters #
####################

print("Launching script with arguments:/n {}".format(args))

print('''
########################
# Load STARsolo output #
########################
''')

# Load adata
adata = sc.read(args.anndata[0])
sample = args.anndata[0].split(".")[0]
print(sample)

print('''
######################
# Filter lowQC cells #
######################
''')

sc.pp.filter_cells(adata, min_genes = 1500) 
sc.pp.filter_genes(adata, min_cells = 3)
adata = adata[adata.obs['percent_mito'] < 0.1, :]

print('''
##################################################
# Remove cell-cycle genes (data-driven approach) #
##################################################
''')

adata.raw = adata.copy()
adata = reptract_utils.per_gene_analysis(adata)

print('''
###############################
# Normalize and log-transform #
###############################
''')

adata = reptract_utils.normalize_log_transform(adata)

print('''
#################################
# HVG selection, PCA, KNN, UMAP #
#################################
''')

adata = reptract_utils.hvgs_pca_umap(adata)

print('''
####################################
# Add scrublet score for each cell #
####################################
''')

#adata = reptract_utils.load_scrublet(adata, sample, args.scrublet[0])

print('''
################################
# Save h5ad file of the sample #
################################
''')

adata.write(args.outdir[0] + sample + ".h5ad")
