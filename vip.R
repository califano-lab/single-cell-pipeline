### takes in a regulon and an expressiond ata set
### computers viper activity, outputs the protein activity matrix

## libraries
library(viper)
## arguments
args <- commandArgs(trailingOnly = TRUE)
datFile <- args[1]
regFile <- args[2]
outFile <- args[3]
## read in data and regulon
expMat <- readRDS(datFile)
reg <- readRDS(regFile)
## compute protein activity and save
pAct <- viper(expMat, regulon = reg, method = 'scale')
saveRDS(pAct, file = outFile)
