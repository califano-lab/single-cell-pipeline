### takes in an expression file and a list of cluster labels
### for each cluster with at least 300 cells, creates a file for that expression set
### stores results in specified directory

## libraries
library(cluster)
## arguments
args <- commandArgs(trailingOnly = TRUE)
datFile <- args[1]
clustLabels <- args[2]
outPath <- args[3]
name <- args[4]
## read in data and cluster labels
dat <- readRDS(datFile)
clustLabels <- readRDS(clustLabels)
## create and save sub expression files
k <- max(clustLabels)
size_thresh <- 300
for (i in 1:k) {
  clust.cells <- which(clustLabels == i)
  if (length(clust.cells) > size_thresh) {
    clust.dat <- dat[,clust.cells]
    save(clust.dat, paste(outPath, name, '_cluster', i, '.rda', sep = ''))
  }
}
