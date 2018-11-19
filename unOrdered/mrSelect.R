### takes in a protein activity matrix
### uses stouffer integration, top/bottom 25 to select 50 MRs for this file
### saves in given output file

## arguments
args <- commandArgs(trailingOnly = TRUE)
pAct.file <- args[1]
outfile <- args[2]
## read in data
pAct <- readRDS(pAct.file)
## stouffer integrate, rank
sInt.proteins <- rowSums(pAct) / sqrt(ncol(pAct))
sInt.proteins <- sort(sInt.proteins)
mrs <- c(sInt.proteins[1:25], tail(sInt.proteins, 25))
## save MRs as table
mrs <- data.frame(mrs); colnames(mrs) <- c('Int_NES')
write.table(mrs, file = outfile, quote = FALSE)
