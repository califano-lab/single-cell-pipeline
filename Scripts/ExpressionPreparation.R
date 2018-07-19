message("==========")
message("Description:")
message("counts/log2(cpm+1) matrix generation")
message("Compulsory Arguments:")
message("--workDir=/full/working/directory/where/results/are/output")
message("--expFile=/full/name/of/raw/counts/file")
message("--jobName=/job/name/")
message("--minCounts=minimun/counts/required")
message("--maxCounts=maximun/counts/required")
message("==========")
message("\n")
#manual

args <- commandArgs(trailingOnly = T)
workDir <- args[grep("--workDir=", args)] 
workDir <- substr(workDir, 11, nchar(workDir))#full working directory
message("working directory:")
message(workDir)
expFile <- args[grep("--expFile=", args)] 
expFile <- substr(expFile, 11, nchar(expFile))#full name of raw expression file
message("raw expression file:")
message(expFile)
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))#job name
message("job name:")
message(jobName)
minCounts <- args[grep("--minCounts=", args)] 
minCounts <- as.numeric(substr(minCounts, 13, nchar(minCounts)))#minimun counts required
message("minimun counts required:")
message(minCounts)
maxCounts <- args[grep("--maxCounts=", args)] 
maxCounts <- as.numeric(substr(maxCounts, 13, nchar(maxCounts)))#maximun counts required
message("maximun counts required:")
message(maxCounts)
#arguments

setwd(workDir)
counts <- get(load(expFile))
counts <- counts[, (colSums(counts) > minCounts & colSums(counts) < maxCounts)]
cpm <- log2(t(t(counts)/(colSums(counts)/1e6)) + 1)
save(cpm, file = paste(jobName, "_cpm.rda", sep = ""))
#data preparation

pdf(paste(jobName, "_Quality_ExpressionPreparation.pdf", sep = ""), width = 12, height = 4)
par(mfrow = c(1, 3))
plot(density(colSums(counts)), main = "# Counts")
plot(density(colSums(counts > 0)), main = "# Genes")
plot(colSums(counts), colSums(counts > 0), xlab = "# Counts", ylab = "# Genes",
     xlim = c(0, 1.2*max(colSums(counts))), ylim = c(0, 1.2*max(colSums(counts > 0))), main = "Saturation")
dev.off()
#QC
message("Finished!")