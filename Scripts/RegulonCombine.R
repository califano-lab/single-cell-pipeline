message("==========")
message("Description:")
message("Generating combined+pruned regulon object (TF, coTF, SIG).")
message("Dependencies:")
message("viper")
message("Compulsory Arguments:")
message("--workDir=/full/working/directory/where/results/are/output")
message("--jobName=/name/of/regulon/file")
message("--netFiles=/full/name/of/included/networks/separated/by/comma")
message("==========")
message("\n")
#manual

args <- commandArgs(trailingOnly = T)
workDir <- args[grep("--workDir=", args)] 
workDir <- substr(workDir, 11, nchar(workDir))#full working directory
message("working directory:")
message(workDir)
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))#name of ARACNe expression file
message("job name:")
message(jobName)
netFiles <- args[grep("--netFiles=", args)]
netFiles <- substr(netFiles, 12, nchar(netFiles))#full name of included networks, seperated by ","
netFiles <- unlist(strsplit(netFiles, split = ","))
message("included networks:")
for (i in 1:length(netFiles)) message(netFiles[i])
#arguments

setwd(workDir)
library(viper)
regul <- lapply(netFiles, function(f){
    message(f)
    regul <- get(load(f))
    regul
})
regul <- do.call(c, regul)
regul <- pruneRegulon(regul)
save(regul, file = paste(jobName, "-all-regulon.rda", sep = ""))
#regulon object
message("Finished!")