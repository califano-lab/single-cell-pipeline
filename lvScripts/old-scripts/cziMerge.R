### cziMerge is a pre-processing step unique to the CZI data set
### this script will take in the acting and resting files, merge them, and save the output matrix
### consider this a 'step 0': the results of this script can the be fed into the true 'step 1', preProcess.R

## libraries
library(stringr)
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-r', '--rest_file'), type="character", help='Input raw count matrix (genesXsamples) for resting samples.'),
  make_option(c('-a', '--act_file'), type="character", help='Input raw count matrix (genesXsamples) for active samples.'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in files and name cells
if(opt$print){ print('Loading data...') }
raw1 <- read.table(opt$rest_file, sep = '\t')
eID <- raw1[,1]; raw1 <- raw1[,-c(1,2)]; rownames(raw1) <- substr(eID, 1, 15)
n1 <- substr(tail(strsplit(opt$rest_file, "/")[[1]], n=1), 1, 5)
colnames(raw1) <- paste(n1, '.', 1:ncol(raw1), sep='')
raw2 <- read.table(opt$act_file, sep = '\t')
eID <- raw2[,1]; raw2 <- raw2[,-c(1,2)]; rownames(raw2) <- substr(eID, 1, 15)
n2 <- substr(tail(strsplit(opt$act_file, "/")[[1]], n=1), 1, 5)
colnames(raw2) <- paste(n2, '.', 1:ncol(raw2), sep='')
## merge into one file, save the raw matrix
if(opt$print){ print('Merging files...') }
shared.genes <- intersect(rownames(raw1), rownames(raw2))
merged.raw <- cbind(raw1[shared.genes,], raw2[shared.genes,])
rm(raw1); rm(raw2); rm(shared.genes)
saveRDS(merged.raw, file=paste(opt$out_dir, opt$out_name, '_mergedRaw.rds', sep=''))