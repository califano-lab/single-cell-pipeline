library(optparse)
## arguments
option_list <- list(
  make_option(c('-i', '--in_file'), type="character", help='Input raw count matrix (genesXsamples) for imputation.'),
  make_option(c('-n', '--out_name'), type="character", help='Output name.'),
  make_option(c('-o', '--out_dir'), type="character", help='Output directory.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
in.dat <- readRDS(opt$in_file)
## make standard distance matrices
euclidean.dist <- dist(t(in.dat), method = 'euclidean')
manhattan.dist <- dist(t(in.dat), method = 'manhattan')
corr.dist <- as.dist(1 - cor(in.dat))
## sims distane matrix
saveRDS(euclidean.dist, file = paste(opt$out_dir, opt$out_name, '_euclidDist.rds', sep = ''))
saveRDS(manhattan.dist, file = paste(opt$out_dir, opt$out_name, '_manhatDist.rds', sep = ''))
saveRDS(corr.dist, file = paste(opt$out_dir, opt$out_name, '_corrDist.rds', sep = ''))