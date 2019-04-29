## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input matrix to be subset (genesXsamples).'),
  make_option(c('-c', '--cluster_labels'), type="character", help='RDS of R cluster object'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
in.dat <- readRDS(opt$input_file)
cluster.labels <- readRDS(opt$cluster_labels)
cluster.labels <- cluster.labels$cluster
## subset the clusters
cluster.table <- table(cluster.labels)
size.thresh <- 300
for (i in 1:length(cluster.table)) {
  if (cluster.table[i] > size.thresh) {
    print(names(cluster.table)[i])
    clust.cells <- names(which(cluster.labels == names(cluster.table)[i]))
    cluster.mat <- in.dat[,clust.cells]
    cluster.mat <- log2(t(t(cluster.mat) / (colSums(cluster.mat) / 1e6)) + 1)
    out.file <- paste(opt$out_dir, opt$out_name, '_cluster-', names(cluster.table)[i], '.rds', sep ='')
    saveRDS(cluster.mat, file = out.file)
    print(dim(cluster.mat))
  }
}
