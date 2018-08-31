message("==========")
message("Description:")
message("UMI matrix generation from 10X Genomics Chromium output")
message("Dependencies:")
message("cellrangerRkit")
message("Compulsory Arguments:")
message("--workDir=/full/working/directory/where/results/are/output")
message("--jobName=/job/name/")
message("==========")
message("\n")
#manual

args <- commandArgs(trailingOnly = T)
workDir <- args[grep("--workDir=", args)] 
workDir <- substr(workDir, 11, nchar(workDir))#full working directory
message("working directory:")
message(workDir)
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))#job name
message("job name:")
message(jobName)
#arguments

library(cellrangerRkit)
setwd(workDir)
umi <- load_cellranger_matrix_from_files("matrix.mtx", "genes.tsv", "barcodes.tsv")
umi <- as.matrix(exprs(umi))
genes <- read.delim("genes.tsv", header=FALSE)
rownames(umi) <- as.character(genes$V2)
umi <- umi[!duplicated(rownames(umi)), ]
save(umi, file = file = paste(jobName, "_umi.rda", sep = ""))
#umi matrix
message("Finished!")
