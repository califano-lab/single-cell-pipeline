### concatenates all protein activity matrices in specified directory into one file
### converts gene names to symbols, subsets for lineage markers
### saves concatenated, filtered file
### requires a specified conversion dictionary

## arguments
args <- commandArgs(trailingOnly = TRUE)
dictFile <- args[1]
pAct.directory <- args[2]
lin.markers <- args[3]
outfile <- args[4]
## read in data
file.names <- dir(pAct.directory, pattern = '*pAct.rds')
pAct <- list()
for (i in 1:length(file.names)) {
  pAct[[i]] <- readRDS(paste(pAct.directory, file.names[i], sep = ''))
}
## concatenate
concat <- pAct[[1]]
for (i in 2:length(pAct)) {
  concat <- rbind(concat, pAct[[i]])
}
## convert ensembl to gene symbol
gene.dict <- read.delim(dictFile, stringsAsFactors = FALSE)
# remove rows with no matches
contained <- which(rownames(concat) %in% gene.dict$Ensembl.ID.supplied.by.Ensembl.)
concat <- concat[contained,]
# replace with gene symbol
row.match <- match(rownames(concat), gene.dict$Ensembl.ID.supplied.by.Ensembl.)
rownames(concat) <- gene.dict[row.match,]$Approved.Symbol
## remove all markers not in the lineage markers
lMarkers <- read.table(lin.markers)
intersect.genes <- intersect(lMarkers[,1], rownames(concat))
concat <- concat[intersect.genes,]
saveRDS(concat, file = outfile)
