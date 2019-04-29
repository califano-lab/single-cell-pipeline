## libraries
library(optparse)
library(stringr)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input matrix to be prepped for scanpy (genesXsamples).'),
  make_option(c('-n', '--out_name'), type="character", help='Name of output.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## load data
input.dat <- readRDS(opt$input_file)
## create annotaiton file and text file
annotation.vect <- substr(colnames(input.dat), 1, 5)
write.table(as.vector(annotation.vect), file = paste(opt$out_dir, opt$out_name, '_scanpy-annotationVect.txt', sep=''),
            sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
write.table(input.dat, file = paste(opt$out_dir, opt$out_name, '_scanpy-inDat.txt', sep = ''), sep = '\t')