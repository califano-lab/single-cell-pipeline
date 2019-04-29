## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-p', '--priority_file'), type="character", help='Input protein activity matrix (proteinsXsamples).'),
  make_option(c('-f', '--fill_file'), type="character", help='Input protein activity matrix (proteinsXsamples).'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
p.mat <- readRDS(opt$priority_file)
f.mat <- readRDS(opt$fill_file)
## copy the fill matrix as the merged matrix
fill.genes <- setdiff(rownames(f.mat), rownames(p.mat))
merged.mat <- rbind(p.mat, f.mat[fill.genes,])
## write out the merged matrix
saveRDS(merged.mat, file = paste(opt$out_dir, opt$out_name, sep=''))
