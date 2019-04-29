#Clustering analysis based on SCANPY pipeline. 
#Implemented by Pasquale Laise

## load libraries
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc
from igraph import Graph

## load arguments
inFile1 = sys.argv[1]
inFile2 = sys.argv[2]
name = sys.argv[3]
outDir = sys.argv[4]

## modify settings 
sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.figdir = outDir
sc.settings.autosave = True
sc.logging.print_versions()
use_first_n_observations = None

## load data
#path = '/Users/pl2659/Documents/SCANPY/'
adata = sc.read(inFile1, cache=False).T  # transpose the data
adata.var_names_make_unique()
#adata
adata.obs['annotation_cell'] = pd.read_csv(inFile2, header=None)[0].values
print(adata.obs)

## normalize data
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

## filter genes, plot dispersion
filter_result = sc.pp.filter_genes_dispersion(
    adata.X, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.filter_genes_dispersion(filter_result, save = '_' + name)

# select filtered genes
adata = adata[:, filter_result.gene_subset]
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)


#do PCA
sc.tl.pca(adata, svd_solver='arpack')
#sc.pl.pca(adata, color='annotation_cell')


sc.pp.neighbors(adata)
sc.logging.print_memory_usage()
sc.tl.umap(adata)
#sc.pl.umap(adata, color='annotation_cell')

sc.tl.louvain(adata,resolution=0.3)
sc.pl.umap(adata, color=['louvain'], save = '_' + name + '_clusters.pdf')
print("Save the figure and proceed with the anlaysis", )
sc.pl.umap(adata, color=['annotation_cell'], save = '_' + name + '_file.pdf')
print("Save the figure and proceed with the anlaysis")
adata.obs.to_csv(outDir + name + '_geneExp-louvainClust.txt', sep='\t')
print("Clustering analysis is done")