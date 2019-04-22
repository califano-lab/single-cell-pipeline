#' PAM clustering over a range of k values.
#' 
#' @param dist.mat A distance matrix for the data to be clustered.
#' @param kmin The minimum k to be used (default of 2).
#' @param kmax The maximum k to be used (default of 5).
#' @param verbose Switch to print status updates
#' @return A list of clustering objects for each value of k.
PamKRange <- function(dist.mat, kmin = 2, kmax = 5, verbose = TRUE) {
  ## packages
  require('cluster')  
  ## generate clusterings for each k
  clusterings <- list()
  for (k in kmin:kmax) {
    if (verbose) { print(paste('Clustering with k=', k, '...', sep = ''))}
    clustering <- pam(dist.mat, k)
    clusterings[[paste('k', k, sep='')]] <- clustering
  }
  return(clusterings)
}

#' Generation of silhouette score for a list of clusters.
#'
#' @param clusterings List of clustering objects.
#' @param dist.mat Distance matrix for the data clustered in clusterings.
#' @param plotPath If specified, will save a plot of the silhouette scores.
#' @return List of silhouette scores.
SilScoreEval <- function(clusterings, dist.mat, plotPath) {
  ## packages
  require('cluster')
  ## find silhouette scores for each cluster
  L <- length(clusterings)
  sil.scores <- c()
  k.vals <- c()
  for (i in 1:L) {
    # identify k for this clustering
    k <- nrow(clusterings[[i]]$medoids); k.vals <- c(k.vals, k)
    # find silhouette score
    sil <- silhouette(clusterings[[i]], dist.mat)
    sil.scores <- c(sil.scores, mean(sil[,3]))
  }
  ## plot if requested
  if (!missing(plotPath)) {
    require('ggplot2')
    plot.dat <- data.frame('k' = k.vals, 'Silhouette.Scores' = sil.scores)
    ggplot(plot.dat, aes(x=k, y=Silhouette.Scores)) + geom_point() + geom_line() +
      ggsave(plotPath, height=2, width=3)
  }
  ## return scores
  return(sil.scores)
}

#' Generates cluster-specific matrices for given data based on a clustering object.
#' 
#' @param dat.mat Data matrix to be split (features X samples).
#' @param clust Clustering object.
#' @param savePath If specified, matrices will be saved. Otherwise, a list of matrices will be returned.
#' @param savePref Preface for file names, if saving.
#' @param sizeThresh Smallest size cluster for which a matrix will be created. Default 300.
#' @return If files are NOT saved, returnes a list of matrices, one for each cluster. Otherwise, returns nothing.
ClusterMatrices <- function(dat.mat, clust, savePath, savePref, sizeThresh = 300) {
  ## set savePath if it is specified
  if(!missing(savePath)) {
    if (!missing(savePref)) {
      savePath <- paste(savePath, savePref, sep = '')
    }
  } else {
    clust.mats <- list()
  }
  ## generate matrices
  clust.table <- table(clust$clustering)
  for (i in 1:length(clust.table)) {
    if (clust.table[i] > sizeThresh) {
      clust.cells <- names(which(clust$clustering == names(clust.table)[i]))
      clust.mat <- dat.mat[, clust.cells]
      if (missing(savePath)) {
        clust.mats[[i]] <- clust.mat
      } else {
        saveRDS(clust.mat, file = paste(savePath, '_', names(clust.table)[i], '.rds', sep = ''))
      }
    }
  }
  ## return if not saving
  if (missing(savePath)) {
    return(clust.mats)
  }
}

#' Generates a pheatmap from the given data and clustering object.
#' 
#' @param dat.mat Data matrix to be used (features X samples).
#' @param clust Clustering object.
#' @param plotTitle Title of the plot.
#' @param plotPath Path to save the plot.
#' @return NULL
ClusterHeatmap <- function(dat.mat, clust, plotTitle, plotPath) {
  require(pheatmap)
  require(RColorBrewer)
  sorted.cells <- sort(clust$cluster)
  pheatmap.mat <- dat.mat[, names(sorted.cells)]
  pheatmap(pheatmap.mat, annotation_col = data.frame('Cluster' = as.factor(clust$cluster)),  
           main = plotTitle, width = 6, height = 8,
           cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(10, 'RdBu')))(100),  fontsize_row = 4,
           filename = plotPath)
}

#' Using a precomputed UMAP, generates a labeled UMAP plot of a clustering.
#' 
#' @param clustering Clustering object.
#' @param umap UMAP object.
#' @param plotTitle Title for the plot.
#' @param plotPath Path to save the plot to.
#' @return NULL
ClusterUMAP <- function(clustering, umap, plotTitle, plotPath) {
  require(ggplot2)
  plot.dat <- data.frame('UMAP1' = umap$layout[,1], 'UMAP2' = umap$layout[,2])
  plot.dat[['cluster']] <- as.factor(clustering$cluster)
  ggplot(plot.dat, aes(x=UMAP1, y=UMAP2, color=cluster)) + geom_point() +
    ggtitle(plotTitle) + ggsave(plotPath, height=8, width=8)
}






