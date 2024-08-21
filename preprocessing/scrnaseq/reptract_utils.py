#!/usr/bin/env python3

import scrublet as scr
import numpy as np
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import os
import sys
import scipy
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42

import seaborn as sns
from scipy import signal
from multiprocessing import Pool
import tables
import scipy.sparse as sp
from typing import Dict

## Reading cellbender outputs can have problem with 3GEX 
# Copy-pasting solution from https://github.com/broadinstitute/CellBender/issues/57#issuecomment-717370745
def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,
                    analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.
    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count matrix.
            True to load a limited set of barcodes: only those analyzed by the
            algorithm. This allows relevant latent variables to be loaded
            properly into adata.obs and adata.obsm, rather than adata.uns.
    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """

    try:
        import anndata
    except ImportError:
        raise ImportError('The anndata package must be installed to use the '
                          'function anndata_from_h5()')

    d = dict_from_h5(file)
    X = sp.csc_matrix((d.pop('data'), d.pop('indices'), d.pop('indptr')),
                      shape=d.pop('shape')).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print('Warning: analyzed_barcodes_only=True, but the key '
                  '"barcodes_analyzed_inds" or "barcode_indices_for_latents" '
                  'is missing from the h5 file. '
                  'Will output all barcodes, and proceed as if '
                  'analyzed_barcodes_only=False')

    # Construct the count matrix.
    adata = anndata.AnnData(X=X,
                            obs={'barcode': d.pop('barcodes').astype(str)},
                            var={'gene_name': (d.pop('gene_names') if 'gene_names' in d.keys()
                                               else d.pop('name')).astype(str)})
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_name', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print('Unable to load data into AnnData: ', key, value, type(value))

    if analyzed_barcodes_only:
        for col in adata.obs.columns[adata.obs.columns.str.startswith('barcodes_analyzed')
                                     | adata.obs.columns.str.startswith('barcode_indices')]:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata

## Function to make the output of anndata_from_h5 the same as sc.read_h5_10x
def uniform_output(adata):
    del adata.obs 
    del adata.uns
    del adata.obsm
    adata.var = adata.var[["id"]]
    adata.var.columns = ["gene_ids"]
    return adata

def read_cellbender_files(filename, raw_file_path, 
               min_n_count = 2000, min_n_gene = 1500, max_n_gene = 8000):
   ## Two-ways of loading cellbender outs
    try:
        adata = sc.read_10x_h5("{r}/{s}.cellbender.out/cellbender_out_filtered.h5".format(r=raw_file_path, s=filename))
    except:
        adata = anndata_from_h5("{r}/{s}.cellbender.out/cellbender_out_filtered.h5".format(r=raw_file_path, s=filename))
        adata = uniform_output(adata)
        
    adata.var_names_make_unique()
    adata.obs_names = [filename+"_"+x.strip("-1") for x in adata.obs_names]
    adata.var["GeneName"] = adata.var_names
    adata.var.columns = ["GeneID", "GeneName"]

    # caculate n_counts / n_genes per cell
    adata.obs['n_counts'] = np.sum(adata.X, axis=1).A1
    adata.obs['n_genes'] = np.sum(adata.X>0,axis=1).A1

    # filter cells
    print("Filtering cells...")
    clist = []
    clist.append(np.array(adata.obs['n_counts'] > min_n_count))
    clist.append(np.array(adata.obs['n_genes'] > min_n_gene))
    clist.append(np.array(adata.obs['n_genes'] < max_n_gene))

    c = np.column_stack(clist).all(axis=1)
    adata = adata[c].copy()

    adata = adata[:,np.argsort(adata.var.GeneID)]
    adata.obs['sample'] = filename
    
    sc.pp.filter_genes(adata, min_cells=3)

    mito_genes = adata.var_names.str.startswith('MT-')
    adata.obs['percent_mito'] = (np.sum(adata.X[:, mito_genes],axis=1).A1) / (np.sum(adata.X,axis=1).A1)
    
    ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL')]
    adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
        
    print("Computing doublets...")
    scrub = scr.Scrublet(adata.X)
    if adata.shape[0] < 30:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=adata.shape[0] - 1)
    else:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets
    #sc.write('/%s/%s'%(outdir, filename),adata)
    return adata

# Function to read starsolo output (when cellbender does not work for a sample) 
def read_starsolo(sample, data_dir,
               min_n_count = 2000, min_n_gene = 1500, max_n_gene = 8000):
    adata = sc.read_10x_mtx(data_dir + sample + '/output/Gene/filtered/', cache=False)
    adata.obs_names = [sample+'_'+i.split('-')[0] for i in adata.obs_names]
    mito_genes = [name for name in adata.var_names if name.startswith('MT-')]
    adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    ribo_genes = [name for name in adata.var_names if name.startswith('RPS') or name.startswith('RPL')]
    adata.obs['percent_ribo'] = np.sum(adata[:, ribo_genes].X, axis = 1).A1 / np.sum(adata.X, axis = 1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    sc.pp.filter_cells(adata, min_genes=min_n_gene)
    sc.pp.filter_cells(adata, max_genes=max_n_gene)
    sc.pp.filter_cells(adata, min_counts=min_n_count)

    print("Computing doublets...")
    scrub = scr.Scrublet(adata.X)
    if adata.shape[0] < 30:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False, n_prin_comps=adata.shape[0] - 1)
    else:
        doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets
    return adata

# Function to plot QC metrics on a per-sample basis 
def qc_plots_sample(adata, sample, outdir):
    plt.tight_layout()
    
    # Set multi-panel figure
    f, axs = plt.subplots(2, 4, figsize=(18, 9))

    # Violin plot of n_genes, n_counts and percent_mito
    v1 = sc.pl.violin(adata, ['n_genes'], show = False, ax=axs[0][0])
    v1.set_title('n_genes violin plot')
    v2 = sc.pl.violin(adata, ['n_counts'], show = False, ax=axs[0][1])
    v2.set_title('n_counts violin plot')
    v3 = sc.pl.violin(adata, ['percent_mito'], show = False, ax=axs[0][2])
    v3.set_title('percent_mito violin plot')

    # Scatterplots of n_counts vs n_genes and n_counts vs percent_mito
    sc.pl.scatter(adata, x='n_counts', y='percent_mito', show = False, title = 'n_counts vs percent_mito scatterplot', ax=axs[0][3])
    sc.pl.scatter(adata, x='n_counts', y='n_genes', show = False, title = 'n_counts vs n_genes scatterplot', ax=axs[1][0])
    
    # Histograms of n_genes and percent_mito
    sns.histplot(adata.obs['n_genes'], bins = 100, kde=True, ax=axs[1][1])
    #plt.axvline(1500, linestyle = '--', color = 'red')
    #plt.axvline(8000, linestyle = '--', color = 'green')

    sns.histplot(adata.obs['percent_mito'], bins = 100, kde=True, ax=axs[1][2])
    #plt.axvline(1500, linestyle = '--', color = 'red')
    #plt.axvline(8000, linestyle = '--', color = 'green')
 
    sns.histplot(adata.obs['percent_ribo'], bins = 100, kde=True, ax=axs[1][3])
    #plt.axvline(1500, linestyle = '--', color = 'red')
    #plt.axvline(8000, linestyle = '--', color = 'green')
    
    f.tight_layout()
    f.savefig(outdir + sample + '_qc_plots.pdf')

# Function to empirically select the minimum number of genes threshold
# NOT CURRENTLY USED 
def choose_min_genes(adata, sample, outdir):
    hist = sns.histplot(adata.obs['n_genes'], bins = 100, kde=True)
    a = hist.get_lines()[0].get_data()[1]

    # maxima : use builtin function to find (max) peaks
    max_peakind = signal.find_peaks_cwt(a, np.arange(1,8))

    # inverse  (in order to find minima)
    inv_data = 1/a
    # minima : use builtin function fo find (min) peaks (use inversed data)
    min_peakind = signal.find_peaks_cwt(inv_data, np.arange(1,8))

    #show results
    print("maxima: ",  max_peakind)
    print("minima: ",  min_peakind)

    my_min_dens = [i for i in min_peakind if i > 10 and i < 30]
    if len(my_min_dens) > 0: 
        my_min_dens = my_min_dens[0]
    else: 
        my_min_dens = min_peakind[1]
    print(my_min_dens)

    my_dict = {'genes' : hist.get_lines()[0].get_data()[0].tolist(), 'density' : hist.get_lines()[0].get_data()[1].tolist()}
    df = pd.DataFrame(my_dict)
    df = df.set_index('density')
    mapping = df['genes'].to_dict()

    my_min_value = mapping[a[my_min_dens]]
    print(my_min_value)

    f, axs = plt.subplots(1, 2, figsize=(8, 4))
    axs[0].plot(a)
    axs[0].axvline(my_min_dens, linestyle = '--', color = 'red')
    axs[0].set_xlabel('KDE values')
    axs[0].set_title('Minimum KDE: ' + str(my_min_dens))

    sns.histplot(adata.obs['n_genes'], bins = 100, kde=True, ax=axs[1])
    plt.axvline(my_min_value, linestyle = '--', color = 'red')
    plt.title('Minimum genes: ' + str(my_min_value).split('.')[0])
    f.tight_layout()
    f.savefig(outdir + sample + '_genes_threshold.pdf')

    return my_min_value

# Function to run gene centric analysis to identify genes that behave like known cell cycle genes 
def per_gene_analysis(adata, species = 'human'):
    bdata = adata.copy()
    # Normalize total counts per cell
    sc.pp.normalize_per_cell(bdata, counts_per_cell_after=1e4)
    # Logarithmize the data matrix
    sc.pp.log1p(bdata)

    # Extract highly variable genes
    sc.pp.highly_variable_genes(bdata)
    highly_variable_genes = bdata.var["highly_variable"]
    bdata = bdata[:, highly_variable_genes]

    # Traspose matrix for a GENE-centered analysis
    bdata = bdata.copy().T

    # Scale data to unit variance and zero mean
    sc.pp.scale(bdata, max_value=10)

    # Scatter plot in PCA coordinates
    sc.tl.pca(bdata)
    bdata.obsm['X_pca'] *= -1  # multiply by -1 to match Seurat
    # Plot the variance ratio
    sc.pl.pca_variance_ratio(bdata, log=True)

    # Compute a neighborhood graph of observations
    sc.pp.neighbors(bdata, n_pcs=20)
    # Embed the neighborhood graph using UMAP
    sc.tl.umap(bdata)
    # Cluster GENES into subgroups using leiden: resolution < 1 to find less clusters
    sc.tl.leiden(bdata, resolution=0.5)

    # Locate ccs cluster
    if species == 'human':
        cycling_genes = ['CDK1','MKI67','CCNB2','PCNA']
        cycling_cdk1 = 'CDK1'
    elif species == 'mouse':
        cycling_genes = ['Cdk1','Mki67']
        cycling_cdk1 = 'Cdk1'

    bdata.obs['known_cyclers'] = [i in cycling_genes for i in bdata.obs_names]
    bdata.obs['known_cyclers'] = bdata.obs['known_cyclers'].astype(int)
    sc.pl.umap(bdata, color=['known_cyclers', 'leiden'], color_map='OrRd')
    print(bdata.obs.loc[[i in cycling_genes for i in bdata.obs_names],'leiden'])
    

    if cycling_cdk1 in bdata.obs_names:
        ccgs_cl = bdata.obs.loc[cycling_cdk1,['leiden']][0]
        print("Cell cycle genes cluster is "+ccgs_cl)
    
         # Add unstructured dict-like annotation for ccgs
        adata.uns['ccgs'] = list(bdata.obs[bdata.obs['leiden']==ccgs_cl].index)
    
        # Remove cc genes
        print('Total number of genes before ccg filter: {:d}'.format(adata.n_vars))
        adata = adata[:,[i not in adata.uns['ccgs'] for i in adata.var_names]]
        print('Total number of genes after ccg filter: {:d}'.format(adata.n_vars))
    else: 
        print("WARNING: CDK1 not present in bdata.obs_names, so not removing cell cycle genes")

    return adata

# Function to normalize and log-transform 
def normalize_log_transform(adata):
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    sc.pp.log1p(adata)
    return adata

# Function to filter HVGs and compute PCA on them 
def hvgs_pca_umap(adata):
    bdata = adata.copy()
    sc.pp.highly_variable_genes(bdata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    for col in ['highly_variable','means', 'dispersions', 'dispersions_norm']:
        adata.var[col] = bdata.var[col]
    bdata = bdata[:, bdata.var['highly_variable']]
    sc.pp.scale(bdata, max_value=10)
    sc.tl.pca(bdata, svd_solver='arpack', n_comps=50)
    #fill NaNs with False so that subsetting to HVGs is possible
    adata.var['highly_variable'].fillna(value=False, inplace=True)
    adata.obsm['X_pca'] = bdata.obsm['X_pca'].copy()
    adata.uns['pca'] = bdata.uns['pca'].copy()
    adata.varm['PCs'] = np.zeros(shape=(adata.n_vars, 50))
    adata.varm['PCs'][adata.var['highly_variable']] = bdata.varm['PCs']
    sc.pl.pca_variance_ratio(adata, log=True)

    # Decide number of PCs to compute neighbourhood graph 
    n_pcs = 25
    sc.pp.neighbors(adata, n_pcs = n_pcs)
    sc.tl.umap(adata)
    return adata

# Benjamini-Hochberg and Bonferroni FDR helper functions.
def bh(pvalues):
    """
    Computes the Benjamini-Hochberg FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    pvalues = np.array(pvalues)
    n = int(pvalues.shape[0])
    new_pvalues = np.empty(n)
    values = [ (pvalue, i) for i, pvalue in enumerate(pvalues) ]
    values.sort()
    values.reverse()
    new_values = []
    for i, vals in enumerate(values):
        rank = n - i
        pvalue, index = vals
        new_values.append((n/rank) * pvalue)
    for i in range(0, int(n)-1):
        if new_values[i] < new_values[i+1]:
            new_values[i+1] = new_values[i]
    for i, vals in enumerate(values):
        pvalue, index = vals
        new_pvalues[index] = new_values[i]
    return new_pvalues

def bonf(pvalues):
    """
    Computes the Bonferroni FDR correction.

    Input:
        * pvals - vector of p-values to correct
    """
    new_pvalues = np.array(pvalues) * len(pvalues)
    new_pvalues[new_pvalues>1] = 1
    return new_pvalues

# Function to run scrublet doublet detection 
def run_scrublet(sample, scrublet_dir, starsolo_dir, file_type):
    #there's loads of clustering going on, so set verbosity low unless you enjoy walls of text
    sc.settings.verbosity = 0  # verbosity: errors (0), warnings (1), info (2), hints (3)

    scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']
    if not os.path.exists(scrublet_dir):
        os.makedirs(scrublet_dir)
    #loop over the subfolders of the rawdata folder
    #import data
    if file_type == 'mtx':
        adata_sample = sc.read_10x_mtx(starsolo_dir + sample + '/',cache=True)
    elif file_type == 'h5':
        adata_sample = sc.read_10x_h5(starsolo_dir + sample + '_filtered_feature_bc_matrix.h5')
    adata_sample.var_names_make_unique()
    #rename cells to SAMPLE_BARCODE
    adata_sample.obs_names = [sample+'_'+i for i in adata_sample.obs_names]
    #do some early filtering to retain meaningful cells for doublet inspection
    sc.pp.filter_cells(adata_sample, min_genes=100)
    sc.pp.filter_genes(adata_sample, min_cells=3)
    #convert to lower to be species agnostic: human mito start with MT-, mouse with mt-
    mito_genes = [name for name in adata_sample.var_names if name.lower().startswith('mt-')]
    # for each cell compute fraction of counts in mito genes vs. all genes
    # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
    adata_sample.obs['percent_mito'] = np.sum(
        adata_sample[:, mito_genes].X, axis=1).A1 / np.sum(adata_sample.X, axis=1).A1
    adata_sample = adata_sample[adata_sample.obs['percent_mito'] < 0.2, :]

    #set up and run Scrublet, seeding for replicability
    np.random.seed(0)
    scrub = scr.Scrublet(adata_sample.X)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(verbose=False)
    adata_sample.obs['scrublet_score'] = doublet_scores

    #overcluster prep. run turbo basic scanpy pipeline
    sc.pp.normalize_per_cell(adata_sample, counts_per_cell_after=1e4)
    sc.pp.log1p(adata_sample)
    sc.pp.highly_variable_genes(adata_sample, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_sample = adata_sample[:, adata_sample.var['highly_variable']]
    sc.pp.scale(adata_sample, max_value=10)
    sc.tl.pca(adata_sample, svd_solver='arpack')
    sc.pp.neighbors(adata_sample)
    #overclustering proper - do basic clustering first, then cluster each cluster
    sc.tl.leiden(adata_sample)
    adata_sample.obs['leiden'] = [str(i) for i in adata_sample.obs['leiden']]
    for clus in np.unique(adata_sample.obs['leiden']):
        adata_sub = adata_sample[adata_sample.obs['leiden']==clus].copy()
        sc.tl.leiden(adata_sub)
        adata_sub.obs['leiden'] = [clus+','+i for i in adata_sub.obs['leiden']]
        adata_sample.obs.loc[adata_sub.obs_names,'leiden'] = adata_sub.obs['leiden']

    #compute the cluster scores - the median of Scrublet scores per overclustered cluster
    for clus in np.unique(adata_sample.obs['leiden']):
        adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_cluster_score'] = \
        np.median(adata_sample.obs.loc[adata_sample.obs['leiden']==clus, 'scrublet_score'])
    #now compute doublet p-values. figure out the median and mad (from above-median values) for the distribution
    med = np.median(adata_sample.obs['scrublet_cluster_score'])
    mask = adata_sample.obs['scrublet_cluster_score']>med
    mad = np.median(adata_sample.obs['scrublet_cluster_score'][mask]-med)
    #let's do a one-sided test. the Bertie write-up does not address this but it makes sense
    zscores = (adata_sample.obs['scrublet_cluster_score'].values - med) / (1.4826 * mad)
    adata_sample.obs['zscore'] = zscores
    pvals = 1-scipy.stats.norm.cdf(zscores)
    adata_sample.obs['bh_pval'] = bh(pvals)
    adata_sample.obs['bonf_pval'] = bonf(pvals)

    #create results data frame for single sample and copy stuff over from the adata object
    scrublet_sample = pd.DataFrame(0, index=adata_sample.obs_names, columns=scorenames)
    for score in scorenames:
        scrublet_sample[score] = adata_sample.obs[score]
    #write out complete sample scores
    scrublet_sample.to_csv(scrublet_dir+sample+'.csv')

# Function to load doublet scores computed with scrublet 
def load_scrublet(adata, sample, scrublet_dir): 
    scorenames = ['scrublet_score','scrublet_cluster_score','zscore','bh_pval','bonf_pval']

    scrdf = pd.read_csv(scrublet_dir + sample + '.csv', header=0, index_col=0)
    scrdf.index = [i.replace('-1', '') for i in scrdf.index]
    for score in scorenames:
        adata.obs[score] = scrdf[score]
    adata.obs['is_doublet'] = adata.obs['bonf_pval'] < 0.01

    # doublets %
    print("Percentage of doublets in dataset: {}".format(adata.obs['is_doublet'].sum() / adata.shape[0]))

    # Fix is_doublet data type to enable plotting correctly 
    adata.obs['is_doublet'] = adata.obs['is_doublet'].astype(int)
    return adata

# Function to make barcode rank plot 
def barcode_plot(sample, starsolo_dir, save_dir):

    # Load unfiltered (raw) data
    adata = sc.read_10x_mtx(starsolo_dir + sample + '/output/Gene/raw/', cache=True)

    # Load filtered data
    fdata = sc.read_10x_mtx(starsolo_dir + sample + '/output/Gene/filtered/', cache=True)

    # Flag cells vs noncells
    adata.obs['is_cell'] = np.where(adata.obs_names.isin(fdata.obs_names.to_list()), 'cell', 'noncell')

    # Compute UMIs per barcode
    adata.obs['umi_num'] = adata.X.sum(axis=1).A

    # Sort the barcodes from higher to lower number of UMIs detected
    df = adata.obs.sort_values(by = ['umi_num'], axis = 0, ascending = False)

    # Add sequential id
    df['barcode_index'] = list(range(0,df.shape[0],1))

    # Make barcode plot and color by cell vs non-cell
    fig, ax = plt.subplots(figsize=(4,4))
    colors = {'cell':'orange', 'noncell':'green'}
    ax.scatter(df['barcode_index'], df['umi_num'], c=df['is_cell'].map(colors))
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('indexed barcodes')
    ax.set_ylabel('UMIs per barcode')
    ax.set_title(sample)
    fig.savefig(save_dir + sample + ".pdf")
