### PACKAGES ###

These are the packages you'll need to run this pipeline:
 - viper
 - cluster
 - ggplot2
 - pheatmap
 - RColorBrewer

### WALKTHROUGH GUIDE ###

This guide will go through a sample analysis using a dummy data set. There are a couple recurring names that will be mentioned:
 - YOUR-CODE-PATH: this is the path to where you have saved the code for this package
 - YOUR-OUT-PATH: this is the path to where you will store results
 - YOUR-OUT-NAME: prefix to be appended to saved files
 - RAW-DAT.rds: the starting file of the analysis, a raw count matrix in .rds format
 - GTEX-NETS.rds: an .rds objects containing the list of GTEx interactomes, if being used.

Start by first loading in the files with the PISCES functions...

```R
source('YOUR-CODE-PATH/process-utils.R')
source('YOUR-CODE-PATH/cluster-functions.R')
source('YOUR-CODE-PATH/viper-utils.R')
```

...as well as the matrix of raw reads:

```R
raw.mat <- readRDS(RAW-DAT.rds)
```

Additionally, we recommend saving intermediate data at each step, since many of these steps will take a decent amount of time with an average sized single cell data set. These savign steps are not including in this walkthrough, but can be achieved with the saveRDS function in R.

### STEP 0: PreProcessing ###

The first step is to generate some simple QC plots based on the raw reads, which will be saved in .pdf format. Next, we filter out low depth cells and low quality genes before converting the raw data into cpm and rank matrices:

```R
QCPlots(raw.mat, 'YOUR-OUT-PATH', 'YOUR-OUT-NAME')
filt.mat <- QCTransform(raw.mat)
cpm.mat <- CPMTransform(filt.mat)
saveRDS(cpm.mat, 'YOUR-OUT-PATH/YOUR-OUT-NAME_cpm.rds')
rank.mat <- RankTransform(cpm.mat
```

By default, the filration step (QCTransform) will remove cells with fewer than 1000 or more than 100000 reads. These numbers are specific to Chromium technology, and should be tailored to your data set.

### STEP 1: R1 Protein Activity and Clustering ###

The first portion of this step is to generate an ARACNe network from the entire data set. This network, refereed to as the *R1-net*, will be used for an intial round of protein activity inference and clustering. 

HOW-TO-GENRATE-NETWORK

With the network created, we will infer protein activity:

```R 
r1.net <- readRDS('YOUR-OUT-PATH/YOUR-OUT-NAME_pruned.rds')
r1.pAct <- viper(rank.mat, r1.net, method = 'none')
```

Optionally, at this step, you can combine the results of the single-cell analysis with those from metaVIPER using the GTEx interactomes. Any regulons not included in the single-cell network will be filled in with data contained in the GTEx interactomes. If you're using this step, do the following:

```R
gtex.nets <- readRDS(GTEX-NETS.rds)
r1.gtex <- viper(rank.mat, gtex.nets, method = 'none')
r1.pAct <- ViperMerge(r1.pAct, r1.gtex)
```

With the first round of protein activity inferred, clustering analysis can be performed on the protein activity-based distance metrix using PAM:

```R
r1.viperDist <- viperSimilarity(r1.pAct)
r1.clust <- PamKRange(r1.viperDist)
r1.clustSil <- SilScoreEval(r1.clust, r1.viperDist)
```

This will generate a set of clusterings for values of *k* between 2 and 5, then generate a vector of average silhouette scores that indicate cluster quality. If you would like the code to automatically generate a plot of the silhouette scores, usign the optional *plotPath* argument with the *SilScoreEval* function. 

For this example, we will assume the optimal clustering occured with *k=2*, but you should analyze your data and choose the clustering with the highest average silhouette score. This clustering will now be subset into matrices that will be used for cluster-specific network generation. Additionally, meta-cells will be inferred from each matrix using the protein activity-based distance matrix.

```R
clust.mats <- ClusterMatrices(cpm.mat, r1.clust$k2)
c1.mCells <- MetaCells(clust.mats[[1]], r1.viperDist)
saveRDS(c1.mCells, file = 'YOUR-CODE-PATH/YOUR-OUT-NAME_r1-c1-mCells.rds')
c2.mCells <- MetaCells(clust.mats[[2]], r1.viperDist)
saveRDS(c2.mCells, file = 'YOUR-CODE-PATH/YOUR-OUT-NAME_r1-c2-mCells.rds')
```

### STEP 2: R2 Network Generation

Before generating networks, we will infer

### STEP 3: R2 Protein Activity and Clustering ###

We will now generate another round of protein activity inference using metaVIPER and the R2 networks:

```R
c1.net <- readRDS('YOUR-OUT-PATH/YOUR-OUT-NAME_c1-r2-net-pruned.rds')
c1.net <- readRDS('YOUR-OUT-PATH/YOUR-OUT-NAME_c2-r2-net-pruned.rds')
r2.pAct <- viper(rank.mat, list('c1' = c1.net, 'c2' = c2.net), method = 'none')
```

Next, clustering is performed using PAM, again using silhouette score to choose the optimal clustering:

```R
r2.viperDist <- viperSimilarity(r2.pAct)
r2.clust <- PamKRange(r2.viperDist)
r2.clustSil <- SilScoreEval(r2.clust, r2.viperDist)
```

### STEP 4: R2 Clustering Analysis ###
























