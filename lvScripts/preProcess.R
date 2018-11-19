### preProcess will prepare raw count data for subsequent analysis
### cells will be filtered based on specified min / max read counts
### genes with no counts will be removed
### filtered matrix will be saved, as will a CPM transformed matrix

## libraries
library(stringr)
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-r', '--raw_file'), type="character", help='Input raw count matrix (genesXsamples).'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.'),
  make_option(c('-a', '--min_count'), type="integer", help='Minimum read count for cells [default]', default=1000),
  make_option(c('-b', '--max_count'), type="integer", help='Maximum read count for cells [default]', default=100000)
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
merged.raw <- readRDS(opt$raw_file)
## filter low quality cells / low count genes, save filtered
if(opt$print){ print('Filtering low quality cells / low count genes...') }
low.thresh <- opt$min_count; high.thresh <- opt$max_count
merged.filtered <- merged.raw[, colSums(merged.raw) > low.thresh & colSums(merged.raw) < high.thresh]
merged.filtered <- merged.filtered[ rowSums(merged.filtered) != 0, ]
rm(merged.raw)
saveRDS(merged.filtered, file=paste(opt$out_dir, opt$out_name, '_mergedFiltered.rds', sep=''))
## save annotation vector
# if(opt$print){ print('Creating annotation vector...') }
# annotation.vect <- substr(colnames(merged.filtered), 1, 5)
# write.table(as.vector(annotation.vect), file = paste(opt$out_dir, opt$out_name, '_annotationVect.txt', sep=''),
#             sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
## convert to CPM, save filtered file
if(opt$print){ print('Converting to CPM...') }
merged.cpm <- log2(t(t(merged.filtered) / (colSums(merged.filtered) / 1e6)) + 1)
saveRDS(merged.cpm, file=paste(opt$out_dir, opt$out_name, '_mergedCPM.rds', sep=''))

