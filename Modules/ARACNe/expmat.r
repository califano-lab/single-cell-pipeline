args <- commandArgs(trailingOnly = T)
wd <- args[1]#working directory
expfn <- args[2]#expression file
acro <- args[3]#acronym
r <- as.numeric(args[4])#minimum expression rate
n <- as.numeric(args[5])#sub-sample
#arguments

expmat <- get(load(expfn))
message(paste("#Samples=", ncol(expmat), sep = ""))
message(paste("N=", n, sep = ""))
if (n > ncol(expmat)){
    message("Using the whole dataset")
}else message(paste("Sub-sample", n, "samples with set.seed", n))
message(paste("Minimum expression rate=", r, sep = ""))
set.seed(n)
expmat <- expmat[, sample(1:ncol(expmat), min(n, ncol(expmat)), replace = F)]
expmat <- expmat[rowSums(expmat > 0) >= r*ncol(expmat), ]
expmat <- cbind(genes=rownames(expmat), round(expmat, 5))
expmat <- rbind(colnames(expmat), expmat)
cat(unlist(apply(expmat, 1, paste, collapse = "\t"), use.names = FALSE), sep = "\n",
    file = paste(wd, "/", acro, "-expmat.dat", sep = ""))
message("Finished!")
#ARACNe expression file
