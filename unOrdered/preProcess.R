### takes in two files; activated and inactivated
### merges, performs CPM, generates a 1000 cell subset for ARACNe
### stores results in specified directory, using the specified name

## libraries
library(stringr)
## arguments
args <- commandArgs(trailingOnly = TRUE)
actFile <- args[1]
restFile <- args[2]
outPath <- args[3]
name <- args[4]
## read in resting; set cell names and prune ensembl names
rest_T <- read.table(restFile, sep = '\t')
ensemble_id <- rest_T[,1]; rest_T <- rest_T[,-c(1,2)]; rownames(rest_T) <- ensemble_id
colnames(rest_T) <- paste('r', name, str_pad(1:ncol(rest_T), 6, pad = '0'), sep = '_')
rownames(rest_T) <- substr(rownames(rest_T), 1, 15)
## read in active; set cell names and prune ensemble names
act_T <- read.table(actFile, sep = '\t')
ensemble_id <- act_T[,1]; act_T <- act_T[,-c(1,2)]; rownames(act_T) <- ensemble_id
colnames(act_T) <- paste('a', name, str_pad(1:ncol(act_T), 6, pad = '0'), sep = '_')
rownames(act_T) <- substr(rownames(act_T), 1, 15)
## merge into one data set
shared_genes <- intersect(rownames(rest_T), rownames(act_T)); length(shared_genes)
merged_raw_counts <- cbind(rest_T[shared_genes,], act_T[shared_genes,])
dim(merged_raw_counts)
## filter for low quality cells / low count genes
low_thresh <- 1000; high_thresh <- 100000
merged_raw_counts.filtered <- merged_raw_counts[, colSums(merged_raw_counts) > low_thresh & colSums(merged_raw_counts) < high_thresh ]
dim(merged_raw_counts.filtered)
merged_raw_counts.filtered <- merged_raw_counts.filtered[ rowSums(merged_raw_counts.filtered) != 0 ,]
dim(merged_raw_counts.filtered)
rm(merged_raw_counts)
saveRDS(merged_raw_counts.filtered, file = paste(outPath, name, '_mergedRaw_filt.rds', sep = ''))
## convert to CPM
merged.cpm <- log2(t(t(merged_raw_counts.filtered) / (colSums(merged_raw_counts.filtered) / 1e6)) + 1)
rm(merged_raw_counts.filtered)
saveRDS(merged.cpm, file = paste(outPath, name, '_mergedCPM.rds', sep = ''))
## create subset for ARACNe
sub_sample <- 1000
names_cells <- sample(colnames(merged.cpm), sub_sample)
expmat <- merged.cpm[, names_cells]
save(expmat, file = paste(outPath, name, '_ARACNeSample.rda', sep = ''))






