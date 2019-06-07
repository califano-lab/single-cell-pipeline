# The pipeline for Protein activity Inference in Single CEllS (PISCES)
---
Authors: Lukas Vlahos, Pasquale Laise, Andrea Califano

Correspondence: Pasquale Laise and Andrea Califano

Contacts: 

* Lukas Vlahos: lv2395@cumc.columbia.edu
* Pasquale Laise: pl2959@cumc.columbia.edu
* Andrea Califano: ac2248@cumc.columbia.edu

---
### Introduction

The pipeline for Protein Activity Inference in Single Cells (PISCES) is a regulatory-network-based methdology for the analysis of single cell gene expression profiles.

PISCES transforms highly variable and noisy single cell gene expression profiles into robust and reproducible protein activity profiles and is centered around two key algorimthms: the Algorithm for the Reconstruction of Accurate Cellular Networks ARACNe [1]; and the algorithm for  Virtual Inference of Protein-activity by Enriched Regulon analysis (VIPER/metaVIPER) [2,3].

Briefly, the ARACNe  algorithm is  one of the most widely used methods for inferring transcriptional interactions from gene expression data. The VIPER algorithm uses the expression of the ARACNe-inferred regulatory targets of a given protein, such as the targets of a transcription factor (TF), as an accurate reporter of its activity. Typically, PISCES  can accurately assess the activity of up to 6000 regulatory proteins  from single cell gene expression profiles,  significantly increasing the ability to analyze the biological function and relevance of gene products whose mRNAs are undetectable in individual cells (e.g. dropout effect).

As currently designed, this pipeline requires a moderate level of computer science know-how. Access to cluster architecture is also advised, as computing networks on a single-core is  time consuming.

### Setup

To run this pipeline, you'll need to have the following packages installed:

* viper
* cluster
* ggplot2
* ggpubr
* umap
* pheatmap
* RColorBrewer
* Matrix
* biomaRt
* psych

```{r include = FALSE}
setwd('C://Users/lvlah/linux/ac_lab/single-cell-pipeline/')
r1.pAct <- readRDS('tutorial/pbmc_r1-pAct.rds')
r2.pAct <- readRDS('tutorial/pbmc_r2-pAct.rds')
```

**NOTE:** This walkthrough will assume that you've set your working directoy to a folder containing the PISCES repo, and that all saved data will be deposited in the same directory. This is not recommended for practical use, and can be changed by specifying full paths when loading or saving in your own applications. 

Start by loading in the PISCES functions as well as the provided test data. Note that this pipeline uses ENSMBL Gene IDs by default:

```{r}
source('functions/process-utils.R')
source('functions/cluster-functions.R')
source('functions/viper-utils.R')
library(ggplot2)
library(ggpubr)
library(viper)
library(pheatmap)
library(RColorBrewer)
raw.mat <- readRDS('tutorial/pbmc.rds')
```

Alternatively, if your data is still in 10x format, you can read it in this manner (you don't need to run this for the purpose of the tutorial):

```{r eval = FALSE}
library(Matrix)
raw.mat <- as.matrix(readMM('data/matrix.mtx'))
genes <- read.table('data/genes.tsv', sep = '\t')
barcodes <- read.table('data/barcodes.tsv', sep = '\t')
colnames(raw.mat) <- barcodes[,1]
rownames(raw.mat) <- genes[,1]
```

We recommend saving intermediate data at each step, since many of these steps will take a considerable amount of time with an average sized single cell data set. These saving steps are not including in this walkthrough, but can be achieved with the saveRDS function in R.

### PreProcessing

First, we perform some cursory QC on the data to check the distribution of read depth and genes detected in each sample:

```{r, fig.align = 'center', fig.width = 9, fig.height = 5}
QCPlots(raw.mat)
```

By default, we remove any genes with no coutns, as well as any samples that have too few (< 1000) or too many (> 100000) UMIs. These thresholds can be adjusted using the arguments of the QCTransform function. Once the data is filtered, we will apply a CPM normalization, then generate a gene expression signature

```{r}
filt.mat <- QCTransform(raw.mat)
cpm.mat <- CPMTransform(filt.mat)
rank.mat <- RankTransform(cpm.mat)
```

By default gene expression signature is generated using a "double rank" approach, which uses the median gene expression of the data set as an internal reference to compute a gene expression signature on a cell by cell basis.

### R1 Network Generation

**NOTE:** Because ARACNe takes a considerable amount of time to run, we recommend setting up a cluster-based implementation. For ease-of-use in this tutorial, we have included the networks generated in this analysis within the tutorial, so that you do not need to generate them yourself. We encourage you to still read an understand this section, but you can also skip directly to **R1 Clustering**.

For a detailed tutorial on how to use ARACNe-AP, consult the tutorial here: https://github.com/califano-lab/ARACNe-AP/blob/master/README.md

By default PISCES generates a single ARACNe network from all the cells in the data set.
This is a default approach that assumes no cells types are known in the data set. If the data set contains multiple cell types and these cell types are known and annotated in the data set (e.g. are experimentally defined  by FACs), you can generate cell-type specific networks. However, we recommend to proceed through the unsupervised approach even if cell types are known ( e.g. experimentally defined), as the unsupervised analysis can either confirm the experimental design and,  potentially, generate novel biological findings. 

The data must first be saved in a format that is compatible with the Java based ARACNe-AP implementation included in this pipeline:

```{r eval = FALSE}
ARACNeTable(cpm.mat, 'tutorial/pbmc-cpm.tsv')
```

ARACNe should be run for each regulator set (TFs, COTFs, and Signaling Proteins). The files containing the lists of the cnadiate master regulator proteins can be found in the */Modules/ARACNe/* directory for both mouse and human data. Once these .tsv's are generated, they should be merged, then combined with the original data to generate a regulon object that can be used to infer protein activity with the VIPER algorithm:

```{r eval = FALSE}
RegProcess('tutorial/pbmc_r1-net-final.tsv', cpm.mat, out.dir = 'tutorial/', out.name = 'pbmc_r1-net-')
```

### R1 Clustering

Once the ARACNe network has been generated, we  can infer protein activity as following:

```{r eval = FALSE}
r1.net <- readRDS('tutorial/pbmc_r1-net-pruned.rds')
r1.pAct <- viper(rank.mat, r1.net, method = 'none')
```

Optionally, at this step, you can combine the results of the single-cell analysis with those from metaVIPER using precomputed GTEx interactomes. 
The rationale here is to further improve the protein activity inference by applying the metaVIPER algorithm. Indeed, metaVPER will make use of GTEX networks to compute the protein activity for the proteins whose regulon (the set of regulatory target genes downstream to a given protein) is not inferred from single cell gene expression profiles.
We suggest loading all the GTEx networks (included in the *GTEx-Nets* directory), but the following exmaple uses only three networks for brevity. Since the GTEx networks utilize Entrez IDs, some gene name conversion will also be necessary.

**NOTE:** As stated, this is an optional (though recommended) step that can be skipped for purposes of this tutorial. Additionally, using all of the avaialble networks is very computationally intensive, and will likely require a cluster-architecture due to the RAM requirements. GTEx networks can be applied only to human samples, this step is not relevant for murine data.

```{r eval = FALSE}
rank.mat.entrez <- Ensemble2Entrez(rank.mat)
lung.net <- get(load('GTEx-Nets/Lung_clean_vst_6cols.rda'))
liver.net <- get(load('GTEx-Nets/Liver_clean_vst_6cols.rda'))
muscle.net <- get(load('GTEx-Nets/Muscle_clean_vst_6cols.rda'))
r1.gtex <- viper(rank.mat.entrez, c(lung.net, liver.net, muscle.net), method = 'none')
r1.gtex <- Entrez2Ensemble(r1.gtex)
r1.pAct <- ViperMerge(r1.pAct, r1.gtex)
```

Once  protein activty has been computed, clustering analysis can be performed by the  PAM algorithm. This pipeline uses *viperSimilarity* as a distance metric and the *PAM* algorithm for clustering, but other methods such as K-Means or louvain can be applied:

```{r, fig.align = 'center'}
r1.viperDist <- as.dist(viperSimilarity(r1.pAct))
r1.clusts <- PamKRange(r1.viperDist, kmin = 2, kmax = 10)
r1.clustSil <- SilScoreEval(r1.clusts, r1.viperDist)
plot.dat <- data.frame('k' = 2:10, 'Silhouette.Scores' = r1.clustSil)
ggplot(plot.dat, aes(x = k, y = Silhouette.Scores)) + geom_point() + geom_line() +
  ggtitle('R1 Clustering Silhouette Scores')
```

This will generate a set of clusters for values of *k* between 2 and 5, then generate a vector of average silhouette scores that indicates cluster quality. It is possible to automatically save a plot of the silhouette scores, usign the optional *plotPath* argument with the *SilScoreEval* function. For the test set used in this tutorial, the optimal number of  clusters is *k=3*. The optimal number of clusters can be defined for each data set looking at the silhouette score. A silhouette score of 0.25 or above is generalyl considered robust [4, 5].



Additionally, in order to improve the quality of the newtorks  by increasing the number of ARACNe regulons, PISCES transforms single cell gene expression profiles into metacell profiles. 
The PISCES meta-cell inference algorithm aims at overcoming the sparse nature of single cell data due to dropouts (inefficient mRNA capture). PISCES accomplishes this step by integrating the expression profiles of cluster of cells (KNN) with similar  protein activty profiles.  After metacell inference, a minimum of two-hundred randomly sampled metacells from each cell type are used to generate cell type specific networks.
All of these steps can be accomplished with the following command (if you are interested in running these steps seperately in your analysis, there are also stand alone commands for each):

```{r eval = FALSE}
r1.clustMats <- MakeCMfA(filt.mat, r1.viperDist, clustering = r1.clusts$k3, out.dir = 'tutorial', out.name = 'pbmc-r1-clusts')
```

### R2 Network Generation

**NOTE:** NOTE: As in the **R1 Network Generation** step, we have  included the ARACNe networks generated in this step within the tutorial. 


The procedure for this step is the same as in **R1 Network Generation**, but must be repeated for each cluster generated in the **R1 Clustering** step. 

```{r eval = FALSE}
c1.net <- RegProcess('tutorial/r2-nets/pbmc-r2-c1_finalNet.tsv', r1.clustMats[[1]], 
                     out.dir = 'tutorial/r2-nets/', out.name = 'pbmc-r2-c1_')
c2.net <- RegProcess('tutorial/r2-nets/pbmc-r2-c2_finalNet.tsv', r1.clustMats[[2]], 
                     out.dir = 'tutorial/r2-nets/', out.name = 'pbmc-r2-c2_')
c3.net <- RegProcess('tutorial/r2-nets/pbmc-r2-c3_finalNet.tsv', r1.clustMats[[3]], 
                     out.dir = 'tutorial/r2-nets/', out.name = 'pbmc-r2-c3_')
```

### R2 Clustering


Cell types-specific  networks  are  then used to infer protein activity  as in the **R1 Clustering** section. Optionally, you can again combine these results with those from the GTEx networks. 

```{r eval = FALSE}
# load in networks
c1.net <- readRDS('tutorial/r2-nets/pbmc-r2-c1_pruned.rds')
c2.net <- readRDS('tutorial/r2-nets/pbmc-r2-c2_pruned.rds')
c3.net <- readRDS('tutorial/r2-nets/pbmc-r2-c3_pruned.rds')
# infer protein activity
r2.pAct <- viper(rank.mat, list('c1' = c1.net, 'c2' = c2.net, 'c3' = c3.net), method = 'none')
```

A new clustering analysis is performed based on protein activty inferred from cell type specific networks

```{r, fig.align = 'center'}
# generate clusterings
r2.viperDist <- as.dist(viperSimilarity(r2.pAct))
r2.clusts <- PamKRange(r2.viperDist, kmin = 2, kmax = 15)
r2.clustSil <- SilScoreEval(r2.clusts, r2.viperDist)
plot.dat <- data.frame('k' = 2:15, 'Silhouette.Scores' = r2.clustSil)
ggplot(plot.dat, aes(x = k, y = Silhouette.Scores)) + geom_point() + geom_line() +
  ggtitle('R2 Clustering Silhouette Scores')
```
In the final step, PISCES  computes a protein acitivity signature of each cluster of cells. In the test set used in this tutorial *k=2* is the optimal clustering, based on silhouette analysi. As stated previously, the the optimal number of clusters  should  defined for each data set data. As a pre-processing step, we convert from ENSG to gene names for ease of interpretation:


```{r}
r2.pAct <- Ensemble2GeneName(r2.pAct)
r2Clust.mrs <- GetMRs(r2.pAct, clustering = r2.clusts$k2$clustering, method = 'Stouffer', numMRs = 50, bottom = FALSE, weights = r2.clusts$k2$silinfo$widths[,3])
```

### Visualization Schemes

We reccomend to visualize the  the differentially activated proteins in each cluster through a heatmap. Indeed, a heatmap can provide a easy visualization of the  differential activty of cluster specific master regulators across all the cells:

```{r, fig.height = 8, fig.width = 6, fig.align = 'center'}
ClusterHeatmap(r2.pAct[ MR_UnWrap(r2Clust.mrs) , ], clust = r2.clusts$k2$clustering, plotTitle = 'R2 Clustering: k = 2')
```

In addition to the heatmap, a uniform manifold approximation and projection (UMAP) can be used  to perform a dimensionality reduction and visualization.
By default, PISCES generate a UMAP using the most representative proteins of each cells. The most representative proteins are defined as the "union" of top active proteins selected on cell by cell basis.

```{r, fig.align = 'center'}
dat.umap <- cbcMR_UMAP(r2.pAct)
ClusterUMAP(r2.clusts$k2$clustering, dat.umap, 'R2 Clustering: k = 2')
```

### Iterative Clustering

In many cases, the first PAM clustering analysis identifies only large clusters of cells which, in turns, can be constitute by subclusters of biological relevance.  
To identify these subcluster if cells we implemented the "IterPAM" function. This functon will perform an iterative clustering analysis untill the silhouette score is higher of given treshold. By default  the treshold on the sihouette score is 0.25, based on the interpretation of  REFs... 

We reccommen to set the parameters of this function according to your data set. If your data has many subpopulations, more iterations may be necessary. If you'd like to allow less robust clusters, the silhouette threshold can be lowered. Different distance functions can also be used, though we recommend the method below if you plan to use protein activity data:

```{r}
dist.func <- function(x) { as.dist(viperSimilarity(x)) }
r2.iClust <- IterPAM(r2.pAct, dist.func)
```

We can visualize the clusters as in  **R2 Clustering**, but with a more detailed annotation about the iterative clustering procedure (annotation of the subclusters):

```{r, fig.height = 8, fig.width = 6, fig.align = 'center'}
iClust.mrs <- GetMRs(r2.pAct, clustering = r2.iClust$clustering[,5], method = 'Stouffer', numMRs = 50, bottom = FALSE)
annot.col <- data.frame(r2.iClust$clustering)
for(i in 1:ncol(annot.col)) { annot.col[,i] <- as.factor(annot.col[,i]) }
pheatmap(r2.pAct[ MR_UnWrap(iClust.mrs) , rownames(annot.col) ], scale = 'row', annotation_col = annot.col,
         color = colorRampPalette(rev(brewer.pal(10, 'RdBu')))(100), height = 8, width = 6,
         cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE, show_rownames = FALSE,
         main = 'R2 Protein Activity Iterative Clustering')
```

A UMAP plot can also be produced in the same manner as above:

```{r, fig.align = 'center'}
ClusterUMAP(r2.iClust$clustering[,5], dat.umap, 'R2 IterClust: i = 5')
```

### Superimposing Markers

The clusters generated by this pipeline are  totally unsupervised. However, the user can  superimpose a set of specific markers for clustering analysis. 
These clusters will provide specific information about the activity of these markers across all the cells.
 
```{r, fig.height = 7, fig.width = 15, fig.align = 'center'}
clust <- r2.iClust$clustering[,2]
plot.dat <- data.frame('UMAP1' = dat.umap$layout[,1], 'UMAP2' = dat.umap$layout[,2])
plot.dat[['cluster']] <- as.factor(clust[rownames(dat.umap$layout)])
clusters <- ggplot(plot.dat, aes(x=UMAP1, y=UMAP2, color=cluster)) + geom_point() + ggtitle('Clusters')
plot.dat[['CCR7']] <- scale(as.numeric(r2.pAct['CCR7',]))
CCR7 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = CCR7)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('CCR7')
plot.dat[['CD8A']] <- scale(as.numeric(r2.pAct['CD8A',]))
CD8A <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = CD8A)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('CD8A')
plot.dat[['MS4A1']] <- scale(as.numeric(r2.pAct['MS4A1',]))
MS4A1 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = MS4A1)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('MS4A1')
plot.dat[['PPBP']] <- scale(as.numeric(r2.pAct['PPBP',]))
PPBP <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = PPBP)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('PPBP')
plot.dat[['MS4A7']] <- scale(as.numeric(r2.pAct['MS4A7',]))
MS4A7 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = MS4A7)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('MS4A7')
plot.dat[['IL7R']] <- scale(as.numeric(r2.pAct['IL7R',]))
IL7R <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = IL7R)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('IL7R')
plot.dat[['CD14']] <- scale(as.numeric(r2.pAct['CD14',]))
CD14 <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = CD14)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('CD14')
plot.dat[['FCER1A']] <- scale(as.numeric(r2.pAct['FCER1A',]))
FCER1A <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = FCER1A)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('FCER1A')
plot.dat[['FCGR3A']] <- scale(as.numeric(r2.pAct['FCGR3A',]))
FCGR3A <- ggplot(plot.dat, aes(UMAP1, UMAP2)) + geom_point(aes(colour = FCGR3A)) + 
  scale_colour_gradientn(colours = c('blue', 'white', 'red')) + ggtitle('FCGR3A')
marker.plot <- ggarrange(clusters, CCR7, CD8A, MS4A1, PPBP, MS4A7, IL7R, CD14, FCER1A, FCGR3A, ncol = 5, nrow = 2)
print(annotate_figure(marker.plot, top = text_grob('R2 Iterclust (i = 2) w/ Canonical PBMC Markers', size = 24)))
```


### References

1.	Lachmann, A., et al., *ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information*. Bioinformatics, 2016. 32(14): p. 2233-5.

2.	Califano, H.D.a.A., *iterClust: Iterative Clustering*. R package version 1.4.0. 2018: https://github.com/hd2326/iterClust.

3.	Ding, H., et al., *Quantitative assessment of protein activity in orphan tissues and single cells using the metaVIPER algorithm*. Nat Commun, 2018. 9(1): p. 1471.

4.  Rosseeuw, P.J., *Journal of Computational and Applied Mathematics* 20 (1987) 53-65

5.  Izenman, A.J., *Modern Multivariate Statistical Techniques. Regression, Classification, and Manifold Learning*. Springer text in statistics, 2008 (Chapter 12)

#### Acknowledgements

Jeremy Dooley - for his advice and expertise in single cell sequencing experiments.

Hongxu Ding - whose work in the Califano laid the groundwork for the development of this pipeline.

Evan Paull - for help with software and tutorial development and testing.
