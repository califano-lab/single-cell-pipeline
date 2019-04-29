## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-a', '--activity_file'), type="character", help='Input protein activity matrix (genesXsamples).'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## read in the data
pAct <- readRDS(opt$activity_file)
## stouffer integrate across the samples
pAct.sInt <- rowSums(pAct) / sqrt(ncol(pAct))
pAct.sInt <- sort(pAct.sInt, decreasing = TRUE)
## get the top and bottom 25 proteins, write to file
top25 <- names(pAct.sInt)[1:25]
bot25 <- tail(names(pAct.sInt), 25)
mrs <- c(top25, bot25)
write.table(mrs, file = paste(opt$out_dir, opt$out_name, '_MRs.txt', sep=''), quote = FALSE, sep = '\t', 
	row.names = FALSE, col.names = FALSE)