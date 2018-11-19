## libraries
library(optparse)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-i', '--input_file'), type="character", help='Input raw count matrix (genesXsamples).'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## load data
if(opt$print){ print('Loading data...') }
raw_mat <- readRDS(opt$input_file)
## generate plots and save
if(opt$print){ print('Generating plots...') }
outfile <- paste(opt$out_dir, opt$out_name, '_QC-Plots.pdf', sep='')
pdf(outfile, onefile = TRUE) 
par(mfrow=c(1,3))
boxplot(colSums(raw_mat), main = "Sequencing depth", frame.plot=F, col="orange")
boxplot(colSums(raw_mat > 0), main = "Detected genes", frame.plot=F, col="cyan")
smoothScatter(colSums(raw_mat), colSums(raw_mat > 0), main = "Saturation plot",
              frame.plot=FALSE, ylab = "Detected genes", xlab = "Sequencing depth", 
              cex = 2, postPlotHook=NULL)
dev.off()
