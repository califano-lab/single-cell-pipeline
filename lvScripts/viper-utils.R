#' Identifies MRs for given data using stouffer integration.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param weights If included, will compute MRs using a weighted stouffer integration.
#' @return Returns the stouffer integrated scores for each protien.
StoufferMRs <- function(dat.mat, weights) {
  # generate dummy weights if missing
  if (missing(weights)) {
    weights = rep(1, ncol(dat.mat))
  }
  # stouffer integrate and return
  sInt <- rowSums(t(t(dat.mat) * weights)) / sqrt(sum(weights ** 2))
  return(sInt)
}

#' Identifies MRs based on ANOVA analysis for a given clustering.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering object.
#' @return A named vector of p-values for each protein
AnovaMRs <- function(dat.mat, clustering) {
  pVals <- c()
  group.vec <- clustering$cluster[colnames(dat.mat)]
  # perform an anova for each protein, storing pValues in a vector
  for (i in 1:nrow(dat.mat)) {
    aov.df <- data.frame('weights' = dat.mat[i,], 'group' = group.vec)
    #print(aov.df)
    aov.test <- aov(weights ~ group, aov.df)
    pVal <- summary(aov.test)[[1]][1,5]
    pVals <- c(pVals, pVal)
  }
  # name and return the vector
  names(pVals) <- rownames(dat.mat)
  return(pVals)
}

#' Returns the master regulators for the given data.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param method 'Stouffer' or 'ANOVA'
#' @param clustering Optional argument for a clustering object. Mandatory for ANOVA method.
#' @param numMRs Number of MRs to return. Default of 50.
#' @param bottom Switch to return downregulated proteins in MR list. Default FALSE>
#' @param weights Optional argument for weights, which can be used in the Stouffer method.
#' @return Returns a list of master regulators, or a list of lists if a clustring is specified.
GetMRs <- function(dat.mat, clustering, method, numMRs = 50, bottom = FALSE, weights, ...) {
  if (method == 'ANOVA') {
    mr.vals <- AnovaMRs(dat.mat, clustering)
  } else if (method == 'Stouffer') {
    # generate dummy weights if not specified
    if (missing(weights)) {
      weights <- rep(1, ncol(dat.mat))
      names(weights) <- colnames(dat.mat)
    }
    # recursive calls for each cluster
    if (missing(clustering)) { # no clustering specified
      mr.vals <- StoufferMRs(dat.mat, weights)
    } else {
      k <- nrow(clustering$medoids)
      mrs <- list()
      for (i in 1:k) {
        # get cluster specific matrix and weights
        clust.cells <- names(which(clustering$cluster == i))
        clust.mat <- dat.mat[, clust.cells]
        clust.weights <- weights[clust.cells]
        # find mrs and add to list
        clust.mrs <- GetMRs(clust.mat, method = method, weights = weights, numMRs = numMRs, bottom = bottom)
        mrs[[paste('c', i, sep = '')]] <- clust.mrs
      }
      return(mrs)
    }
  } else {
    print('Invalid method: must be "Stouffer" or "ANOVA".')
  }
  # return appropriate portion of MR list
  mr.vals <- sort(mr.vals, decreasing = TRUE)
  if (bottom) {
    return(c(mr.vals[1:numMRs], tail(mr.vals, numMRs)))
  } else {
    return(mr.vals[1:numMRs])
  }
}

#' Generates a UMAP based on a set of proteins.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param mrs List of proteins to use as master regulators.
#' @return UMAP object.
MRUMAP <- function(dat.mat, mrs) {
  require(umap)
  dat.mat <- dat.mat[mrs,]
  dat.umap <- umap(t(dat.mat))
  return(dat.umap)
}

#' Generates a meta cell matrix for given data.
#' 
#' @param dat.mat Raw gene expression matrix (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor inference.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 5.
#' @param subSize If specified, number of metaCells to be subset from the final matrix. No subsetting occurs if not incldued.
#' @return A matrix of meta cells (genes X samples).
MetaCells <- function(dat.mat, dist.mat, numNeighbors = 5, subSize) {
  # prune distance matrix if necessary
  dist.mat <- dist.mat[colnames(dat.mat), colnames(dat.mat)]
  # KNN function
  KNN <- function(dist.mat, numNeighbors){
    dist.mat <- as.matrix(dist.mat)
    n <- nrow(dist.mat)
    neighbor.mat <- matrix(0L, nrow = n, ncol = k)
    for (i in 1:n) {
      neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
    }
    return(neighbor.mat)
  }
  knn.neighbors <- KNN(dist.mat, 10)
  # create imputed matrix
  imp.mat <- matrix(0, nrow = nrow(dat.mat), ncol = ncol(dat.mat))
  rownames(imp.mat) <- rownames(dat.mat); colnames(imp.mat) <- colnames(dat.mat)
  for (i in 1:ncol(dat.mat)) {
    neighbor.mat <- dat.mat[,c(i, knn.neighbors[i,])]
    imp.mat[,i] <- rowSums(neighbor.mat)
  }
  if (missing(subSize)) {
    return(imp.mat)
  } else {
    return(imp.mat[, sample(colnames(imp.mat), min(subSize, ncol(imp.mat)))])
  }
}

#' Merges two viper matrices, giving priority to one over the other.
#' 
#' @param p.mat Priority viper matrix (proteins X samples). Proteins here will override those in the other matrix.
#' @param q.mat Secondary viper matrix (proteins X samples). Proteins here will fin in for gaps in the priority matrix.
#' @return A merged viper matrix.
ViperMerge <- function(p.mat, q.mat) {
  fill.genes <- setdiff(rownames(q.mat), rownames(p.mat))
  merged.mat <- rbind(p.mat, q.mat[fill.genes,])
  return(merged.mat)
}






