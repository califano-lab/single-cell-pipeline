### takes in a data matrix (protein activity or gene expression)
### generates hierarchical clustering from k=2 to 10
### for each clustering, generates a PCA plot, silhouette plot, and a pheatmap plot
### stores clustering and plots in specified directory, preficed with given name

## libraries
library(ggplot2)
library(pheatmap)
library(cluster)
library(RColorBrewer)
## arguments
args <- commandArgs(trailingOnly = TRUE)
datFile <- args[1]
outPath <- args[2]
name <- args[3]
## read in data, generate PCA and distance matrix
dat <- readRDS(datFile)
dat.pca <- prcomp(t(dat))
saveRDS(dat.pca, file = paste(outPath, name, '-datPCA.rds', sep = ''))
dat.dist <- dist(t(dat), method = 'manhattan')
saveRDS(dat.dist, file = paste(outPath, name, '-datDist.rds', sep = ''))
## generate clustering
hClust <- hclust(dat.dist, method = 'ward')
saveRDS(hClust, file = paste(outPath, name, '-hClust.rds', sep = ''))
## generate plots for each k
silScores <- c()
for (k in 2:10) {
  clust <- cutree(hClust, k = k)
  saveRDS(clust, file = paste(outPath, name, '-hClust-k', k, '.rds', sep = ''))
  plotTitle <- paste("HClust of ", name, ', k=', k, sep="")
  ## set colors
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cluster.colors <- gg_color_hue(k); names(cluster.colors) <- 1:k
  ## pca plot
  plot.dframe <- as.data.frame(dat.pca$x)
  plot.dframe$'Clusters' <- as.factor(clust)
  ggplot(plot.dframe, aes(PC1, PC2, colour = Clusters)) +
    geom_point() + 
    ggtitle(plotTitle)
  ggsave(paste(outPath, name, '-pca-hClust-k', k, '.jpg', sep = ''), plot = last_plot(), height = 7, width = 7, unit = 'in')
  ## silhouette plot
  clust.sil <- silhouette(clust, dat.dist)
  silScores <- c(silScores, mean(clust.sil[,3]))
  jpeg(paste(outPath, name, '-sil-hClust-k', k, '.jpg', sep = ''), width = 700, height = 700)
  plot(clust.sil, col = cluster.colors, border = NA, main = plotTitle)
  dev.off()
  ## pheatmap
  sorted.labels <- data.frame('cluster' = sort(clust))
  pheatmap(dat[,rownames(sorted.labels)], cluster_cols = FALSE, cluster_rows = TRUE, labels_col = '', show_rownames = TRUE,
           annotation_col = sorted.labels, annotation_colors = list('cluster' = cluster.colors),
           main = plotTitle, fontsize_row = 5, height=6, width=5, color = colorRampPalette(rev(brewer.pal(10, 'RdBu')))(100),
           filename = paste(outPath, name, '-pheatmap-hClust-k', k, '.jpg', sep = ''))
}
## create a plot of silhouette scores
jpeg(paste(outPath, name, '-silScores', '.jpg', sep = ''), width = 700, height = 700)
plot(2:10, silScores, xlab = 'k', ylab = 'Silhouette Score', 
     main = paste('Silhouette Scores of ', name, sep = ''))
lines(2:10, silScores)
dev.off()








