## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input matrix to be subset (genesXsamples).'),
  make_option(c('-c', '--cluster_labels'), type="character", help='Text file with cluster and sample labels.'),
  make_option(c('-s', '--sample_index'), type="integer", help='Index of sample indices in cluster label file.'),
  make_option(c('-l', '--cluster_index'), type="integer", help='Index of cluster labels in cluster lable file.'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
in.dat <- readRDS(opt$input_file)
cluster.labels <- read.table(opt$cluster_labels, header = TRUE, stringsAsFactors = FALSE)
s.ind <- opt$sample_index
l.ind <- opt$cluster_index
## subset the clusters
cluster.table <- table(cluster.labels[,l.ind])
size.thresh <- 300
for (i in 1:length(cluster.table)) {
  if (cluster.table[i] > size.thresh) {
    print(names(cluster.table)[i])
    ind.match <- which(cluster.labels[,l.ind] == names(cluster.table)[i])
    cluster.samples <- cluster.labels[ind.match, s.ind]
    cluster.mat <- in.dat[,cluster.samples]
    out.file <- paste(opt$out_dir, opt$out_name, '_cluster-', names(cluster.table)[i], '.rds', sep ='')
    saveRDS(cluster.mat, file = out.file)
    print(dim(cluster.mat))
  }
}










