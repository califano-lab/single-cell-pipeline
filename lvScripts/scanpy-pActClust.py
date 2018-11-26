## libraries
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc
from igraph import Graph
## arguments
inFile1 = sys.argv[1]
inFile2 = sys.argv[2]
name = sys.argv[3]
outDir = sys.argv[4]
## adjust figure save settings
sc.settings.figdir = outDir
sc.settings.autosave = True
## load data
adata = sc.read(inFile1, cache=False).T 
adata.obs['annotation_cell'] = pd.read_csv(inFile2, header=None)[0].values
## filter based on dispersion
filter_result = sc.pp.filter_genes_dispersion(
    adata.X, flavor='seurat', n_top_genes=100, log=False)
adata = adata[:, filter_result.gene_subset]
## plot PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='annotation_cell', save = '_' + name + '_file.pdf')
## plot UMAP of data
sc.pp.neighbors(adata, metric="correlation")
sc.logging.print_memory_usage()
sc.tl.umap(adata)
sc.pl.umap(adata, color='annotation_cell', save = '_' + name + '_file.pdf')
## plot UMAP of clusters and CD2
sc.tl.louvain(adata,resolution=0.3)
sc.pl.umap(adata, color=['louvain'], save = '_' + name + '_clusters.pdf')
#sc.pl.umap(adata, color=['CD2'], save = '_' + name + '_CD2.pdf')
## write the UMAP for downstream analysis in R
np.savetxt(outDir + name + '_pAct-louvainClust.txt', adata.obsm.X_umap, delimiter='\t')