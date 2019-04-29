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

#' Iterative clustering using PAM
#' 
#' @param dat.mat Matrix of data (features X samples).
#' @param dist.func Functiont to compute a distance matrix from the data.
#' @param iter.max Maximum number of iterations to use
#' @param sil.thresh Minimum cluster quality threshold (default of 0.25).
#' @param iter.num Internal parameter for tracking iteration depth.
#' @return A data frame with samples in rows and iterations in columns, containing cluster labels for each iteration.
IterPAM <- function(dat.mat, dist.func, iter.max = 3, sil.thresh = 0.25, iter.num = 1) {
  print(iter.num)
  # generate distance matrix
  dist.mat <- dist.func(dat.mat)
  # perform clustering
  iter.clust <- PamKRange(dist.mat, verbose = FALSE)
  sil.scores <- unlist(lapply(iter.clust, function(x) {x$silinfo$avg.width} ))
  best.sil <- which.max(sil.scores)
  # identify best clustering, prepare return object
  opt.clust <- iter.clust[[names(best.sil)]]
  clust.df <- data.frame(opt.clust$clustering); colnames(clust.df) <- c(paste0('iter.', iter.num))
  opt.k <- length(table(opt.clust$clustering))
  # return NULL if the clustering does not pass the threshold (ie is homogeneous)
  if (max(sil.scores) < sil.thresh) { return(NULL) }
  # return if this is the final iteration
  if (iter.num >= iter.max) { return(clust.df) }
  # iterate to next layer of clustering
  clust.iters <- list()
  for (i in 1:opt.k) {
    # cluster the samples for this cluster
    clust.samps <- which(opt.clust$clustering == i)
    sub.clust <- IterPAM(dat.mat[,clust.samps], dist.func, 
                         iter.max = iter.max, iter.num = iter.num + 1, sil.thresh = sil.thresh)
    clust.iters[[i]] <- sub.clust
  }
  # if all the sub clusterings were null, return the df as is
  if (all(lapply(clust.iters, is.null))) { return(clust.df) }
  # build data frame with sub clusters
  iter.cols <- max(unlist(lapply(clust.iters, function(x) { ncol(x) })))
  while (ncol(clust.df) < iter.cols + 1) { # add dummy columns, that will be filled in later
    clust.df <- cbind(clust.df, clust.df[,1])
  }
  space.vect <- rep(0, iter.cols)
  for (i in 1:length(clust.iters)) { # add in the sub clustering labels
    if (!is.null(clust.iters[[i]])) { # for all sub clusterings that were not null
      # build to correct size
      while(ncol(clust.iters[[i]]) < iter.cols) {
        clust.iters[[i]] <- cbind(clust.iters[[i]], clust.iters[[i]][,1])
      }
      for (j in 1:iter.cols) {
        clust.df[rownames(clust.iters[[i]]), (1 + j)] <- clust.iters[[i]][, j] + space.vect[[j]]
        space.vect[[j]] <- space.vect[[j]] + length(unique(clust.iters[[i]][,j]))
      }
    } else {
      space.vect <- space.vect + 1
    }
  }
  # rename and return
  colnames(clust.df) <- paste0('iter.', iter.num + 0:(ncol(clust.df) - 1))
  return(clust.df)
}
