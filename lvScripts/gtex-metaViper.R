### gtex-metaViper will run metaViper using the GTEx bulk ARACNe networks
### this ahs a dedicated script (for now) because the GTEx networks use Entrez and our data uses ENSMBL

## libraries
library(stringr)
library(optparse)
library(viper)
## arguments
option_list <- list(
  make_option(c('-q', '--quiet'), action='store_false', default=TRUE, help='Suppresses status updates.', dest='print'),
  make_option(c('-c', '--cpm_file'), type="character", help='Input cpm matrix (genesXsamples).'),
  make_option(c('-i', '--net_dir'), type="character", help='Directory with the interactomes to use.'),
  make_option(c('-k', '--convert_dict'), type="character", help='Dictionary for converting between gene names.'),
  make_option(c('-n', '--out_name'), type="character", help='Output files name convention.'),
  make_option(c('-d', '--out_dir'), type="character", help='Directory for output.')
)
opt <- parse_args(OptionParser(option_list = option_list))
## load data and comput double rank signature
if(opt$print){ print('Loading data and calculating signature...') }
cpm.data <- readRDS(opt$cpm_file)
rank.data <- apply(cpm.data, 2, rank)
median <- apply(rank.data, 1, median)
mad <- apply(rank.data, 1, mad)
rank.signature <- (rank.data - median) / mad
rm(cpm.data); rm(rank.data)
## convert to entrez ids
if(opt$print){ print('Converting gene names...') }
gene.dict <- read.delim(opt$convert_dict, header = TRUE, stringsAsFactors = FALSE)
contained.ensembl <- rownames(rank.signature)[which(rownames(rank.signature) %in% gene.dict$Ensembl.Gene.ID)]
subDict <- gene.dict[which(gene.dict$Ensembl.Gene.ID %in% contained.ensembl),]
subDict <- subDict[which(subDict$Entrez.Gene.ID != ''),]
rank.signature <- rank.signature[subDict$Ensembl.Gene.ID,]
rownames(rank.signature) <- subDict$Entrez.Gene.ID
## load regulons in a list
if(opt$print){ print('Loading networks...') }
file_names=as.list(dir(path = opt$net_dir, pattern="*.rda"))
regList <- list()
for(i in 1:length(file_names)) 
{
  new.reg <- get(load(paste0(opt$net_dir,file_names[[i]]))) 
  new.name <-paste0("reg_",i)
  regList[[new.name]] <- new.reg
}
## perform metaviper 
if(opt$print){ print('Performing metaVIPER...') }
vip.mat <- viper(rank.signature, regulon = regList, method = 'none')
## convert back to ENSEMBL
rownames(vip.mat) <- subDict[match(rownames(vip.mat), subDict$Entrez.Gene.ID),]$Ensembl.Gene.ID
## save results
saveRDS(vip.mat, file=paste(opt$out_dir, opt$out_name, '_GTExActivity.rds', sep=''))