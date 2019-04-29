library(cluster)
library(ggplot2)
library(optparse)
library(umap)
library(pheatmap) 
library(RColorBrewer)
library(biomaRt)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input matrix to be clustered (proteinsXsamples).'),
  make_option(c('-m', '--kmin'), type='integer', default=2, help='Minimum k to test, default of 2.'),
  make_option(c('-k', '--kmax'), type='integer', default=5, help='Maximum k to test, default of 5'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
dat.mat <- readRDS(opt$input_file)
kmin <- opt$kmin; kmax <- opt$kmax
out_name <- opt$out_name; out_dir <- opt$out_dir

dat.mat <- readRDS('C://Users/lvlah/linux/ac_lab/data/czi/czi-r2/d1-boneMarrow/d1-BM_cluster-net-viper.rds')
kmin <- 2; kmax <- 5
out_name <- 'd1-boneMarrow_r2-reClust'
out_dir <- 'C://Users/lvlah/linux/ac_lab/data/czi/reClust/d1-boneMarrow/'

## fucntion to convert gene names
Ensemble2GeneName<-function(dataset2Convert) {
  require(biomaRt)
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  names_dataset<-getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'),filters = 'ensembl_gene_id', values = (rownames(dataset2Convert)), mart = ensembl)
  rownames(names_dataset)<-make.unique(names_dataset$ensembl_gene_id)
  dataset2Convert_2<-merge(names_dataset,dataset2Convert,by=c("row.names"))
  dim(dataset2Convert_2)
  rownames(dataset2Convert_2)<-make.unique(dataset2Convert_2$hgnc_symbol)
  GeneName_dataset<-dataset2Convert_2[,-c(1:4)]
  head(GeneName_dataset[,1:4])
  return(GeneName_dataset)
}
## generate UMAP for the data
set.seed(1)
dat.umap <- umap(t(dat.mat))
plot.dat <- data.frame('UMAP1' = dat.umap$layout[,1], 'UMAP2' = dat.umap$layout[,2])
saveRDS(dat.umap$layout, file = paste(out_dir, out_name, '_r2-UMAP.rds', sep = ''))
## generate distance matrix for the data
library(viper)
dat.dist <- as.dist(viperSimilarity(dat.mat))
saveRDS(dat.dist, file = paste(out_dir, out_name, '_r2-vipSim-dist.rds', sep = ''))
## generate clusterings for the range of k's specified, producing silhouette plots as well
allScores <- c()
mrActivity <- c()
for (k in kmin:kmax) {
  ## clustering and silhouette score updating
  print(k)
  clustering <- pam(dat.dist, k)
  sil <- silhouette(clustering$cluster, dat.dist)
  allScores <- c(allScores, mean(sil[,3]))
  ## umap plot
  plot.dat[['cluster']] <- as.factor(clustering$cluster)
  ggplot(plot.dat, aes(x=UMAP1, y=UMAP2, color=cluster)) + geom_point() +
    ggtitle(paste(out_name, ': k=', k, sep = '')) +
    ggsave(paste(out_dir, out_name, '_k', k, '-UMAP.jpeg', sep = ''), height=8, width=8)
  ## silhouette plot
  jpeg(paste(out_dir, out_name, '_k', k, '-silhouette.jpeg', sep = ''))
  plot(sil, border = NA)
  dev.off()
  ## save the clustering
  saveRDS(clustering, file = paste(out_dir, out_name, '_k', k, '-clustering.rds', sep=''))
  ## mr_analysis (stouffer integration)
  mrs <- data.frame(row.names = rownames(dat.mat))
  avgMR <- 0
  sig <- c()
  mr.thresh <- 50; nP <- 50
  for (i in 1:k) {
    # subset the cluster matrix and generate the stouffer integration
    clust.cells <- names(which(clustering$cluster == i))
    clust.mat <- dat.mat[,clust.cells]
    pAct.sInt <- rowSums(clust.mat) / sqrt(ncol(clust.mat))
    # store the stouffer integration in the mrs data frame, identify the significant mrs for cluster selection
    mrs <- cbind(mrs, pAct.sInt)
    avgMR <- avgMR + (mean(pAct.sInt[pAct.sInt > mr.thresh]) * length(clust.cells))
    # signature for the pheatmap
    pAct.sInt <- sort(pAct.sInt, decreasing = TRUE)
    sig <- c(sig, pAct.sInt[1:nP])
  }
  colnames(mrs) <- 1:k
  write.csv(mrs, file = paste(out_dir, out_name, '_k', k, '_masterRegulators.csv', sep = ''), quote=FALSE)
  mrActivity <- c(mrActivity, avgMR / ncol(dat.mat))6
  ## pheatmap
  sig <- unique(sig)
  sorted.cells <- sort(clustering$cluster)
  pheatmap.mat <- dat.mat[sig, names(sorted.cells)]; pheatmap.mat <- Ensemble2GeneName(pheatmap.mat)
  pheatmap(pheatmap.mat, annotation_col = data.frame('Cluster' = as.factor(clustering$cluster)),  
           main = paste(out_name, ': k=', k, sep = ''), width = 6, height = 8,
           cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE,
           color = colorRampPalette(rev(brewer.pal(10, 'RdBu')))(100),  fontsize_row = 4,
           filename = paste(out_dir, out_name, '_k', k, '-signature-heatmap.jpeg', sep= ''))
}
## generate plot of the silhouette scores
plot.dat2 <- data.frame('k' = kmin:kmax, 'Silhouette.Scores' = allScores, 'Signature.Activity' = mrActivity)
ggplot(plot.dat2, aes(x=k, y=Silhouette.Scores)) + geom_point() + geom_line() +
  ggtitle(paste(out_name, ': Silhouette Scores', sep = '')) + 
  ggsave(paste(out_dir, out_name, '_silhouette-scores.jpeg', sep = ''), height=8, width=8)
## generate plot of average MR activity
ggplot(plot.dat2, aes(x=k, y=Signature.Activity)) + geom_point() + geom_line() +
  ggtitle(paste(out_name, ': Signature MR Average Activity', sep = '')) + 
  ggsave(paste(out_dir, out_name, '_signatue-activity.jpeg', sep = ''), height=8, width=8)

