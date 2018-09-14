#Step 4 of single cell pipeline
#SampleLevelActivityAnalysis.r
#input file - *_cpm.rda, *.txt (annoFile.txt)
#output file - *_exp.rda,*_pca_exp.rda, *_mds_pca_exp.rda, *_Quality_SampleLevelActivityAnalysis.pdf, *_Sample_SampleLevelactivityAnalysis.pdf, *_Anno_SampleLevelActivityAnalysis.pdf, *_iterClust_SampleLevelActivityAnalysis.pdf

#####--------------------------------------Section 1-------------user messages------------------------------------------########
rmarkdown::render("SampleLevelActivityAnalysis.R", "pdf_document")
message("==========")
message("Description:")
message("Activity analysis oat sample level,
including PCA/MDS dimension reduction analysis, QC report generation,
cell type annotation and iterClust analysis.")
message("Dependencies:")
message("viper, cluster, iterClust, made4")#libraries
message("Compulsory Arguments:")
message("--workDir=/full/working/directory")
message("--sampleFiles=/full/name/of/included/samples/separated/by/comma")
message("--netFiles=/full/name/of/included/networks/separated/by/comma")
message("--jobName=/job/name/")
message("Optional Arguments:")
message("--annoFile=/full/name/of/expression/markers/file")#optional
message("tab-delimited .txt with cell names in 1st column and gene names in 2nd column, see sample file")
message("==========")
message("\n")
#manual
#####---------------------------------End of--Section 1----------user messages------------------------------------------########

#####--------------------------------------Section 2-------------Sample files-------------------------------------------########
#arguments
args <- commandArgs(trailingOnly = T)
workDir <- args[grep("--workDir=", args)] 
workDir <- substr(workDir, 11, nchar(workDir))#full working directory
message("working directory:")
message(workDir)
#files
sampleFiles <- args[grep("--sampleFiles=", args)]
sampleFiles <- substr(sampleFiles, 15, nchar(sampleFiles))#full name of included samples, seperated by ","
sampleFiles <- unlist(strsplit(sampleFiles, split = ","))
message("included samples:")
#networks
for (i in 1:length(sampleFiles)) message(sampleFiles[i])
netFiles <- args[grep("--netFiles=", args)]
netFiles <- substr(netFiles, 12, nchar(netFiles))#full name of included networks, seperated by ","
netFiles <- unlist(strsplit(netFiles, split = ","))
message("included networks:")
for (i in 1:length(netFiles)) message(netFiles[i])
jobName <- args[grep("--jobName=", args)] 
jobName <- substr(jobName, 11, nchar(jobName))#job name
message("job name:")
message(jobName)
if(length(grep("--annoFile=", args)) != 0){#optional argument, annotation
    annoFile <- args[grep("--annoFile=", args)] 
    annoFile <- substr(annoFile, 12, nchar(annoFile))#full name of expression marker file
}else annoFile <- NA
message("expression marker file:")
message(annoFile)
#arguments
#####----------------------------------End of ----Section 2-------Sample files------------------------------------------########

#####--------------------------------------Section 3-------------Data Preparation---------------------------------------########
#viper/metaviper
setwd(workDir)
library(viper)
source("/ifs/scratch/c2b2/ac_lab/st3179/sc-rnaseq-pipeline/test-demo/single-cell-pipeline-master/Modules/metaVIPER/metaVIPER.R")
source("/ifs/scratch/c2b2/ac_lab/st3179/sc-rnaseq-pipeline/test-demo/single-cell-pipeline-master/Modules/ColorGradient/ColorGradient.R")

gene <- unique(unlist(lapply(sampleFiles, function(f){
    exp <- get(load(f))
    rownames(exp)
})))
exp <- lapply(sampleFiles, function(f, gene){
    message(f)
    exp <- get(load(f))
    f <- sapply(strsplit(f, split = "/"), function(x) x[length(x)])
    f <- unlist(strsplit(f, split = "_cpm.rda"))
    t <- matrix(0, nrow(exp), ncol(exp), dimnames = list(gene, paste(f, 1:ncol(exp), sep = "_")))
    t[rownames(exp), ] <- exp
    t
}, gene=sort(gene, decreasing = T))#sort
exp <- do.call(cbind, exp)
rank <- apply(exp, 2, rank)
median <- apply(rank, 1, median)
mad <- apply(rank, 1, mad)
rank <- (rank - median)/mad
regulon <- lapply(netFiles, function(f){#regulon
    message(f)
    regul <- get(load(f))
    regul
})#ranking
#metaviper
vp <- metaVIPER(eset = rank, regulon = regulon, weight = "mean", method = "none")
save(vp, file = paste(jobName, "_vp.rda", sep = ""))#output
pca_vp <- prcomp(t(vp))
save(pca_vp, file = paste(jobName, "_pca_vp.rda", sep = ""))
mds_vp <- cmdscale(as.dist(viperSimilarity(vp)))
save(mds_vp, file = paste(jobName, "_mds_vp.rda", sep = ""))
#data preparation
#####----------------------------------End of----Section 3---------Data Preparation---------------------------------------########

#####--------------------------------------Section 4-------------Quality control------------------------------------------########
library(made4)
cat <- getcol(21)[c(-1, -3, -12)]
pdf(paste(jobName, "_Quality_SampleLevelActivityAnalysis.pdf", sep = ""), width = 12, height = 6)
par(mfrow = c(2, 4))
#mds
plot(mds_vp, xlab = "MDS-1", ylab = "MDS-2", cex = 0.5, pch = 16,
     col = cat[as.factor(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1]))], main = "Samples")
legend("topright", legend = paste(levels(as.factor(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1]))),
                                  table(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1])), sep = ":"),
       fill = cat, bty = "n", border = NA)
plot(density(colSums(exp > 0)), main = "# Genes")
plot(mds_vp, xlab = "MDS-1", ylab = "MDS-2", col = geneColor(exp), cex = 0.5, pch = 16, main = "# Genes")
barplot(signif(pca_vp$sdev[1:10]^2/sum(pca_vp$sdev^2)*100, 2), xlab = "PC", ylab = "%Variance", main = "")
#pca
for (i in 2:5){
    plot(pca_vp$x[, 1], pca_vp$x[, i], col = cat[as.factor(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1]))],
         xlab = paste("PC1, ", signif(pca_vp$sdev[1]^2/sum(pca_vp$sdev^2)*100, 2), "%", sep = ""),
         ylab = paste("PC", i, ", ", signif(pca_vp$sdev[i]^2/sum(pca_vp$sdev^2)*100, 2), "%", sep = ""), cex = 0.5, pch = 16, main = "")
}
dev.off()
#data visualization & QC

pdf(paste(jobName, "_Sample_SampleLevelActivityAnalysis.pdf", sep = ""), width = 12,
    height = (floor(length(table(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1])))/4)+1)*3)
par(mfrow = c(floor(length(table(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1])))/4)+1, 4))
for (x in levels(as.factor(sapply(strsplit(rownames(mds_vp), split = "_"), function(x) x[1])))){
    smoothScatter(mds_vp[grep(x, rownames(mds_vp)), ], nrpoints = 0, nbin = 256,
                  xlim = range(mds_vp[, 1]), ylim = range(mds_vp[, 2]), xlab = "MDS-1", ylab = "MDS-2", main = x)
}
dev.off()
#individual datasets
#####---------------------------------End of--Section 4-----------Quality control------------------------------------------########

#####--------------------------------------Section 5-------------Handling annotation---------------------------------------########
if (!is.na(annoFile)){
    anno <- read.delim(annoFile, header=FALSE)
    #annoFile
    pdf(paste(jobName, "_Anno_SampleLevelActivityAnalysis.pdf", sep = ""), width = 16, height = (floor(2*nrow(anno)/4)+1)*4)
    par(mfrow = c(floor(2*nrow(anno)/4)+1, 4))
    #mds
    for (i in 1:nrow(anno)){
        plot(mds_vp, xlab = "MDS-1", ylab = "MDS-2", cex = 0.5, pch = 16,
             col = expColor(anno[i, 2], exp, c("Black", "Red")), main = paste("Exp.", anno[i, 1], anno[i, 2]))
        legend("bottomleft", legend = paste("exp=", signif(range(exp[anno[i, 1], ]), 2), sep = ""),
               fill = c("Black", "Red"), bty = "n", border = NA)
        plot(mds_vp, xlab = "MDS-1", ylab = "MDS-2", cex = 0.5, pch = 16,
             col = vpColor(anno[i, 2], vp), main = paste("Act.", anno[i, 1], anno[i, 2]))
        legend("bottomleft", legend = c("Act.<-5", "Act.=0", "Act.>5"), fill = c("Blue", "Grey", "Red"), bty = "n", border = NA)
    }
    dev.off()
    #figure
}
#annotation
#####----------------------------------End of----Section 5--------Handling annotation---------------------------------------########

#####--------------------------------------Section 6-------------Iterclust--------------------------------------------------########
#clustering
library(iterClust)
library(cluster)
library(made4)
library(viper)
cat <- getcol(21)[c(-1, -3, -12)]
CC <- function (dset, iteration){
    dist <- as.dist(viperSimilarity(dset))
    range <- seq(2, 5, by = 1)
    clust <- vector("list", length(range))
    for (i in 1:length(range)) clust[[i]] <- pam(dist, range[i])$clustering
    return(clust)
}#range
CE <- function (dset, iteration, clust){
    dist <- as.dist(viperSimilarity(dset))
    clustEval <- vector("numeric", length(clust))
    for (i in 1:length(clust)) {
        clustEval[i] <- mean(silhouette(clust[[i]], dist)[, "sil_width"])
    }
    return(clustEval)
}#evaluation
OE <- function (dset, clust, iteration){
    dist <- as.dist(viperSimilarity(dset))
    obsEval <- vector("numeric", length(clust))
    return(silhouette(clust, dist)[, "sil_width"])
}
CH <- function (clustEval, iteration) return(clustEval > 0 * iteration + 0.01)
OO <- function (obsEval, iteration) return(obsEval < 0 * iteration)
FS <- function (dset, iteration, feature) return(rownames(dset))
clust_vp <- iterClust(dset = structure(vp, dimnames=list(rownames(vp), 1:ncol(vp))),
                      maxIter = 4, coreClust = CC, clustEval = CE, featureSelect = FS,
                      obsEval = OE, clustHetero = CH, obsOutlier = OO)
save(clust_vp, file = paste(jobName, "_clust_vp.rda", sep = ""))
pdf(paste(jobName, "_iterClust_SampleLevelActivityAnalysis.pdf", sep = ""), width = 10, height = 10)
par(mfrow = c(2, 2))
for (j in 1:4){
    COL <- structure(rep(1, ncol(vp)), names = 1:ncol(vp))
    for (i in 1:length(clust_vp$cluster[[j]])) COL[clust_vp$cluster[[j]][[i]]] <- cat[i]
    plot(mds_vp, col = COL, cex = 0.5, pch = 16, xlab = "MDS-1", ylab = "MDS-2", main = paste("iterClust, Round", j))
    legend("topright", legend = c(paste("C", 1:i, sep = ""), "Outlier"), fill = c(cat[1:i], 1), bty = "n", border = NA)
}
dev.off()
#iterClust
message("Finished!")
#####-----------------------------------End of---Section 6-------------Iterclust--------------------------------------------------########
