#Step 1 of single cell pipeline
#ExpressionPreparation.r
#input file - *.rda/umi.rda, raw count matrix file
#output file - *cpm.rda and quality control results (*_Quality_ExpressionPreparation.pdf)

#####--------------------------------------Section 1-------------user messages------------------------------------------########
#rmarkdown::render("ExpressionPreparation.R", "pdf_document")
#help for end user
message("==========")
message("Description:")
message("counts/log2(cpm+1) matrix generation")
message("Compulsory Arguments:") #compulsory args
message("--workDir=/full/working/directory/where/results/are/output")
message("--expFile=/full/name/of/raw/counts/file")
message("--jobName=/job/name/")
message("--minCounts=minimun/counts/required") #cell filtering parameter
message("--maxCounts=maximun/counts/required")
message("==========")
message("\n")
#manual
#####---------------------------------End of--Section 1----------user messages------------------------------------------########

#####--------------------------------------Section 2-------------arguments----------------------------------------------########
#starting arguments
args <- commandArgs(trailingOnly = T)
workDir <- args[grep("--workDir=", args)] #working directory 
workDir <- substr(workDir, 11, nchar(workDir))#full working directory
#general messages
message("working directory:")
message(workDir)
#expression files
expFile <- args[grep("--expFile=", args)]#arguments for expression files 
expFile <- substr(expFile, 11, nchar(expFile))#full name of raw expression file
message("raw expression file:")
message(expFile)
#job specification
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))#job name
message("job name:")
message(jobName)
#count calculation
minCounts <- args[grep("--minCounts=", args)]#minimum count calculation 
minCounts <- as.numeric(substr(minCounts, 13, nchar(minCounts)))#minimun counts required
message("minimun counts required:")
message(minCounts)
maxCounts <- args[grep("--maxCounts=", args)]#maximum count calulcation 
maxCounts <- as.numeric(substr(maxCounts, 13, nchar(maxCounts)))#maximun counts required
message("maximun counts required:")
message(maxCounts)
#arguments
#####---------------------------------End of-----Section 2-------arguments----------------------------------------------########

#####--------------------------------------Section 3-------------Calculation of cpm-------------------------------------########
#cpm matrix creation
setwd(workDir)
#working with counts
counts <- get(load(expFile))#expression file loading
counts <- counts[, (colSums(counts) > minCounts & colSums(counts) < maxCounts)]
cpm <- log2(t(t(counts)/(colSums(counts)/1e6)) + 1)#cpm calculation
save(cpm, file = paste(jobName, "_cpm.rda", sep = ""))
#data preparation

#####---------------------------------End of-----Section 3-------Calulation of cpm--------------------------------------########
#####--------------------------------------Section 4-------------output charts------------------------------------------########
#output
#density chanrts
pdf(paste(jobName, "_Quality_ExpressionPreparation.pdf", sep = ""), width = 12, height = 4)
par(mfrow = c(1, 3))
plot(density(colSums(counts)), main = "# Counts")#counts plot
plot(density(colSums(counts > 0)), main = "# Genes")#genes plot
plot(colSums(counts), colSums(counts > 0), xlab = "# Counts", ylab = "# Genes",
     xlim = c(0, 1.2*max(colSums(counts))), ylim = c(0, 1.2*max(colSums(counts > 0))), main = "Saturation")#saturation plot
dev.off()
#QC
message("Finished!")
#####---------------------------------End of-----Section 4------output charts---------------------------------------------########
