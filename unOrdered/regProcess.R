### takes in the path of an expresion file and a network
### converts the regulon using aracne2regulon, making it usable for Viper
### saves it in the specified file

## libraries
library(viper)
## arguments
args <- commandArgs(trailingOnly = TRUE)
afile <- args[1]
print(afile)
eset.file <- args[2]
print(eset.file)
output <- args[3]
## read in data and convert
eset <- get(load(eset.file))
processed.reg <- aracne2regulon(afile = afile, eset = eset, format = '3col')
## prune if size is specified
if (length(args) == 4) {
  processed.reg <- pruneRegulon(processed.reg, as.integer(args[4]), adaptive = FALSE, eliminate = TRUE)
}
## output
saveRDS(processed.reg, file = output)
