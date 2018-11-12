#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#"""
#Created on Mon Oct 29 13:54:58 2018

#@author: pl2659
#"""
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
import scanpy.api as sc
from igraph import Graph

inFile1 = sys.argv[1]
inFile2 = sys.argv[2]
outFile= sys.argv[3]

adata = sc.read(inFile1, cache=False).T 

adata.obs['annotation_cell'] = pd.read_csv(inFile2, header=None)[0].values


#i f you do not normalize the data, the software will automatiaccly use the protein acity data 


filter_result = sc.pp.filter_genes_dispersion(
    adata.X, flavor='seurat', n_top_genes=100, log=False)

adata = adata[:, filter_result.gene_subset]

#plot PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='annotation_cell')

sc.pp.neighbors(adata, metric="correlation")
sc.logging.print_memory_usage()
sc.tl.umap(adata)
sc.pl.umap(adata, color='annotation_cell')

sc.tl.louvain(adata,resolution=0.3)
sc.pl.umap(adata, color=['louvain'])
sc.pl.umap(adata, color=['CD2'])

#write the UMAP for downstream analysis in R
np.savetxt(outFile, adata.obsm.X_umap, delimiter='\t')

