#' Identifies MRs for given data using stouffer integration.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param cluster Vector of cluster lables. If not included, integrates the entire matrix.
#' @param weights A named vector of sample weights. If included, stouffer integration is weighted.
#' @return Returns the stouffer integrated scores for each protien.
StoufferMRs <- function(dat.mat, cluster, weights) {
  # generate dummy weights if missing
  if (missing(weights)) {
    weights <- as.numeric(rep(1, ncol(dat.mat))); names(weights) <- colnames(dat.mat)
  }
  # perform integration across full matrix if cluster was missing
  if (missing(cluster)) {
    sInt <- rowSums(t(t(dat.mat) * weights))
    sInt <- rowSums(t(t(dat.mat) * weights)) / sqrt(sum(weights ** 2))
    return(sort(sInt, decreasing = TRUE))
  }
  # separate cluster specific matrices
  k <- length(table(cluster))
  mrs <- list()
  for (i in 1:k) { # for each cluster
    clust.cells <- names(cluster)[which(cluster == i)]
    clust.mat <- dat.mat[, clust.cells]
    clust.weights <- weights[clust.cells]
    clust.mrs <- StoufferMRs(clust.mat, weights = clust.weights)
    mrs[[paste('c', i, sep = '')]] <- sort(clust.mrs, decreasing = TRUE)
  }
  return(mrs)
}

#' Identifies MRs based on ANOVA analysis for a given clustering.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Clustering vector
#' @return A named vector of p-values for each protein
AnovaMRs <- function(dat.mat, clustering) {
  pVals <- c()
  group.vec <- clustering[colnames(dat.mat)]
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

#' Performs a bootstrap t-test between two sample vectors x and y. Returns a log p-value.
#' 
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogBootstrapTTest <- function(x, y, bootstrap.num = 100) {
  x.n <- length(x); y.n <- length(y)
  log.pValue <- c()
  ## perform test for each bootstrap
  for (i in 1:bootstrap.num) {
    # create bootstraps
    x.boot <- sample(1:x.n, size = x.n, replace = TRUE)
    x.boot <- x[x.boot]
    y.boot <- sample(1:y.n, size = y.n, replace = TRUE)
    y.boot <- y[y.boot]
    # perform t.test
    test.res <- t.test(x = x.boot, y = y.boot, alternative = "two.sided")
    # generate log p-value
    log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
    log.pValue <- c(log.pValue, log.p)
  }
  # return mean log p-value
  return(mean(log.pValue))
}

#' Performs a t-test between two sample vectors x and y. Returns a log p-value.
#' @param x Vector of test values.
#' @param y Vector of reference values.
#' @param bootstrap.num Number of bootstraps to use. Default of 100.
#' @return A signed log p-value.
LogTTest <- function(x, y) {
  test.res <- t.test(x, y, alternative = 'two.sided')
  log.p <- 2*pt(q = abs(test.res$statistic), df = floor(test.res$parameter), log.p = TRUE, lower.tail = FALSE)*(-sign(test.res$statistic))
  return(log.p)
}

#' Identifies MRs based on a bootstraped Ttest between clusters.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param clustering Vector of cluster labels.
#' @param bootstrap.num Number of bootstraps to use. Default of 10 
#' @return Returns a list of lists; each list is a vector of sorted log p-values for each cluster.
BTTestMRs <- function(dat.mat, clustering, bootstrap.num = 100) {
  # set initial variables
  clustering <- clustering
  k <- length(table(clustering))
  mrs <- list()
  # identify MRs for each cluster
  for (i in 1:k) {
    print(paste('Identifying MRs for cluster ', i, '...', sep = ''))
    mrs.mat <- matrix(0L, nrow = nrow(dat.mat), ncol = bootstrap.num)
    rownames(mrs.mat) <- rownames(dat.mat)
    # split test and ref matrices
    clust <- names(table(clustering))[i]
    clust.vect <- which(clustering == clust)
    test.mat <- dat.mat[, clust.vect]; ref.mat <- dat.mat[, -clust.vect]
    t.n <- ncol(test.mat); r.n <- ncol(ref.mat)
    # for each bootstrap
    for (b in 1:bootstrap.num) {
      test.boot <- test.mat[, sample(colnames(test.mat), size = t.n, replace = TRUE)]
      ref.boot <- ref.mat[, sample(colnames(ref.mat), size = t.n, replace = TRUE)]
      # for each gene
      for (g in rownames(dat.mat)) {
        mrs.mat[g, b] <- LogTTest(test.boot[g,], ref.boot[g,])
      }
    }
    # sort and add to list
    mList <- sort(rowMeans(mrs.mat), decreasing = TRUE)
    mrs[[clust]] <- mList
  }
  # return 
  return(mrs)
}

#' Returns the master regulators for the given data.
#' 
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param method 'Stouffer' or 'ANOVA'
#' @param clustering Optional argument for a vector of cluster labels.
#' @param numMRs Number of MRs to return per cluster. Default of 50.
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
      k <- length(table(clustering))
      mrs <- list()
      for (i in 1:k) {
        # get cluster specific matrix and weights
        clust.cells <- names(which(clustering == i))
        clust.mat <- dat.mat[, clust.cells]
        print(dim(clust.mat))
        clust.weights <- weights[clust.cells]
        # find mrs and add to list
        clust.mrs <- GetMRs(clust.mat, method = method, weights = clust.weights, numMRs = numMRs, bottom = bottom)
        print(head(clust.mrs))
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

#' Identifies MRs on a cell-by-cell basis and returns a merged, unique list of all such MRs.
#'
#' @param dat.mat Matrix of protein activity (proteins X samples).
#' @param numMRs Default number of MRs to identify in each cell. Default of 25.
#' @return Returns a list of master regulators, the unique, merged set from all cells.
CBCMRs <- function(dat.mat, numMRs = 25) {
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:numMRs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # return
  return(cbc.mrs)
}

#' Make Cluster Metacells for ARACNe. Will take a clustering and produce saved meta cell matrices.
#'
#' @param dat.mat Matrix of filtered, un-normalzied gene expression (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor calculation. We recommend using a viper similarity matrix.
#' @param numNeighbors Number of neighbors to use for each meta cell; if 0, meta cellsare not computed.. Default of 5.
#' @param clustering Vector of cluster labels. 
#' @param subSize Size to subset the data too. Since 200 cells is adequate for ARACNe runs, this allows for speedup. Default of 200.
#' @param pseudo If TRUE, uses a pseudo count for normalization. Default of FALSE.
#' @param out.dir Directory for sub matrices to be saved in.
#' @param out.name Optional argument for preface of file names. 
MakeSubMats <- function(dat.mat, dist.mat, numNeighbors = 0, clustering, subSize = 200, pseudo = FALSE, out.dir, out.name = '') {
  # generate cluster matrices
  clust.mats <- ClusterMatrices(dat.mat, clustering)
  # produce sub matrices and save for each cluster matrix
  k <- length(clust.mats)
  sub.mats <- list()
  for (i in 1:k) { # for each cluster
    mat <- clust.mats[[i]]
    if (numNeighbors > 0) { # metacells if specified
      mat <- MetaCells(mat, dist.mat, numNeighbors, subSize)
    } else {
      mat <- mat[, sample(colnames(mat), subSize)]
    }
    # cpm transformation
    mat <- CPMTransform(mat, pseudo = pseudo)
    sub.mats[[i]] <- mat
    # save file
    file.name <- paste(out.dir, out.name, '_clust-', i, '.rds', sep = '')
    saveRDS(mat, file = file.name)
  }
  return(sub.mats)
}

#' Generates a meta cell matrix for given data.
#' 
#' @param dat.mat Raw gene expression matrix (genes X samples).
#' @param dist.mat Distance matrix to be used for neighbor inference.
#' @param numNeighbors Number of neighbors to use for each meta cell. Default of 10.
#' @param subSize If specified, number of metaCells to be subset from the final matrix. No subsetting occurs if not incldued.
#' @param cpm Switch to control transformation to cpm. Default TRUE.
#' @return A matrix of meta cells (genes X samples).
MetaCells <- function(dat.mat, dist.mat, numNeighbors = 10, subSize, cpm = TRUE) {
  # prune distance matrix if necessary
  dist.mat <- as.matrix(dist.mat)
  dist.mat <- dist.mat[colnames(dat.mat), colnames(dat.mat)]
  dist.mat <- as.dist(dist.mat)
  # KNN function
  KNN <- function(dist.mat, k){
    dist.mat <- as.matrix(dist.mat)
    n <- nrow(dist.mat)
    neighbor.mat <- matrix(0L, nrow = n, ncol = k)
    for (i in 1:n) {
      neighbor.mat[i,] <- order(dist.mat[i,])[2:(k + 1)]
    }
    return(neighbor.mat)
  }
  knn.neighbors <- KNN(dist.mat, numNeighbors)
  # create imputed matrix
  imp.mat <- matrix(0, nrow = nrow(dat.mat), ncol = ncol(dat.mat))
  rownames(imp.mat) <- rownames(dat.mat); colnames(imp.mat) <- colnames(dat.mat)
  for (i in 1:ncol(dat.mat)) {
    neighbor.mat <- dat.mat[,c(i, knn.neighbors[i,])]
    imp.mat[,i] <- rowSums(neighbor.mat)
  }
  # subset if specified
  if (!missing(subSize)) {
    s.size <- min(subSize, ncol(imp.mat))
    imp.mat <- imp.mat[, sample(colnames(imp.mat, s.size)) ]
  }
  # cpm if specified
  if (cpm) {
    imp.mat <- CPMTransform(imp.mat)
  }
  # return
  return(imp.mat)
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

#' Processes ARACNe results into a regulon object compatible with VIPER.
#'
#' @param a.file ARACNe final network .tsv.
#' @param exp.mat Matrix of expression from which the network was generated (genes X samples).
#' @param out.dir Output directory for networks to be saved to.
#' @param out.name Optional argument for prefix of the file name.
RegProcess <- function(a.file, exp.mat, out.dir, out.name = '.') {
  require(viper)
  processed.reg <- aracne2regulon(afile = a.file, eset = exp.mat, format = '3col')
  saveRDS(processed.reg, file = paste(out.dir, out.name, 'unPruned.rds', sep = ''))
  pruned.reg <- pruneRegulon(processed.reg, 50, adaptive = FALSE, eliminate = TRUE)
  saveRDS(pruned.reg, file = paste(out.dir, out.name, 'pruned.rds', sep = ''))
}

#' Performs a stouffer integration of a list of viper matrices.
#' 
#' @param vip.mats List of viper matrices (proteins X samples).
#' @param weights Vector of weights for the stouffer integration. If not included, all matrices are weighted equally.
#' @return An integrated viper matrix (proteins X samples).
VIPIntegrate <- function(vip.mats, weights) {
  # set weights, if not specified
  if (missing(weights)) {
    weights <- rep(1, length(vip.mats))
  }
  # identify set of regulators and create final matrix
  regs <- unique(Reduce(union, lapply(vip.mats, rownames)))
  int.mat <- matrix(0L, nrow = length(regs), ncol = ncol(vip.mats[[1]]))
  rownames(int.mat) <- regs; colnames(int.mat) <- colnames(vip.mats[[1]])
  # integrate each reg
  for (i in 1:length(regs)) {
    # generate matrix of all instances of this regulon
    reg.mat <- as.data.frame( lapply(vip.mats, function(x) {
      if (regs[i] %in% rownames(x)) {
        return( x[ regs[i] ,] )
      } else {
        return( rep(NA, ncol(x)) )
      }
    } ) )
    incl.vec <- which(!is.na(reg.mat[1,]))
    # integrata
    c.weights <- weights[incl.vec]
    sInt <- rowSums( t(t(reg.mat[, incl.vec]) * c.weights) ) / sqrt(sum(c.weights ** 2))
    int.mat[regs[i] ,] <- sInt
  }
  # return
  return(int.mat)
}

#' MetaVIPER implementation that will perform a weighted stouffer integration based on highest NES.
#' 
#' @param ges Gene Expression Signature (features X samples)
#' @param net.list List object with the networks to be used
#' @param use.nets Optional argument to sslect the top n networks. If not specified, all networks are used.  
#' @param ret.weights Optional argument to return the network weight matrix as well as the VIPER matrix. FALSE by default.
#' @return Either a viper matrix, or a list with a viper matrix and the network weight matrix.
WeightedVIPER <- function(ges, net.list, use.nets, ret.weights = FALSE) {
  num.nets <- length(net.list)
  num.samps <- ncol(ges)
  ## create weight matrix
  w.mat <- matrix(0L, nrow = num.nets, ncol = ncol(ges))
  colnames(w.mat) <- colnames(ges); rownames(w.mat) <- names(net.list)
  ## run VIPER with each network
  print('Generating VIPER matrices...')
  vip.list <- list()
  for (i in 1:num.nets) {
    vip.list[[i]] <- viper(ges, net.list[i], method = 'none')
  }
  names(vip.list) <- names(net.list)
  ## count for each gene
  print('Generating weights...')
  uni.genes <- unique(unlist(lapply(vip.list, rownames)))
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(0)
        }}))
      max.ind <- which.max(abs(nes.vals))
      w.mat[max.ind, s] <- w.mat[max.ind, s] + 1
    }
  }
  ## integration
  print('Integrating...')
  int.mat <- matrix(0L, nrow = length(uni.genes), ncol = num.samps)
  rownames(int.mat) <- uni.genes; colnames(int.mat) <- colnames(ges)
  for (g in uni.genes) {
    for (s in 1:num.samps) {
      nes.vals <- unlist(lapply(vip.list, function(x){ 
        if (g %in% rownames(x)) {
          return(x[g,s])
        } else {
          return(NA)
        }}))
      w.vals <- w.mat[,s][!is.na(nes.vals)]
      w.vals <- w.vals / sum(w.vals)
      nes.vals <- nes.vals[!is.na(nes.vals)]
      # if use.nets are specified, subset to the top n (or use all if n > length)
      if (!missing(use.nets)) {
        w.order <- order(w.vals, decreasing = TRUE)
        w.vals <- w.vals[ w.order[1:min(length(w.order), use.nets)] ]
        nes.vals <- nes.vals[ w.order[1:min(length(nes.vals), use.nets)] ]
      }
      int.mat[g,s] <- sum(nes.vals * w.vals) / sqrt(sum(w.vals**2))
    }
  }
  ## return
  if (ret.weights) {
    return( list('viper' = int.mat, 'weights' = w.mat) )
  } else {
    return( int.mat )
  }
}