## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input matrix to be subset (genesXsamples).'),
  make_option(c('-c', '--convert_dict'), type="character", help='Dictionary for gene name conversion.'),
  make_option(c('-s', '--start_index'), type="integer", help='Index of sample indices in cluster label file.'),
  make_option(c('-l', '--dest_index'), type="integer", help='Index of cluster labels in cluster lable file.'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in data
in.dat <- readRDS(opt$input_file)
convert.dict <- readRDS(opt$convert_dict)
## convert names
if(opt$print){ print('Converting gene names...') }
contained.names <- rownames(in.dat)[which(rownames(in.dat) %in% convert.dict[,opt$start_index])]
subDict <- convert.dict[which(convert.dict[,opt$start_index] %in% contained.names),]
subDict <- subDict[which(subDict[,opt$dest_index] != ''),]
in.dat <- in.dat[as.character(subDict[,opt$start_index]),]
rownames(in.dat) <- subDict[,opt$dest_index]
## save matrix with converted names
saveRDS(in.dat, file = paste(opt$out_dir, opt$out_name, sep = ''))
l