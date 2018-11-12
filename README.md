

## A single cell pipeline for the identification of cell type-specific master regulators
####### by Pasquale Laise (last update 10-31-2018)

Paths for the HPC for datasets
* Datasets: 
  * ```/ifs/scratch/c2b2/ac_lab/CZI/peter_data```
  
* ARACNe repository: 
  * ```/ifs/scratch/c2b2/ac_lab/CZI/aracne-repository```

## Example: DONOR1 (Lung)

````
rm(list=ls(all=T))
library(cluster)
library(data.table)

````
### Load data sets

Load "resting T cells" matrix

````
rest_T<-fread("~/PP001swap.filtered.matrix.txt.bz2",sep="\t",colClasses = "numeric"))
rest_T<-as.data.frame(rest_T)
class(rest_T[2,3])
ensemble_id<-rest_T[,1]
colnames(rest_T)<-paste("PP001",seq(1:length(colnames(rest_T))),sep=".")
dim(rest_T)
str(rest_T)
head(rest_T[,1:4])
rest_T<-rest_T[,-c(1,2)]
rownames(rest_T)<-ensemble_id
head(rest_T[,1:4])
annot_rest<-rep("PP001",length(colnames(rest_T)))

````
Load  "activated T cells" matrix

````
act_T<-fread("~/PP002swap.filtered.matrix.txt.bz2",sep="\t",colClasses = "numeric")
act_T<-as.data.frame(act_T)
dim(act_T)
str(act_T)
ensemble_id<-act_T[,1]
colnames(act_T)<-paste("PP002",seq(1:length(colnames(act_T))),sep=".")
act_T<-act_T[,-c(1,2)]
rownames(act_T)<-ensemble_id
head(act_T[,1:4])
annot_act_T<-rep("PP002",length(colnames(act_T)))
dim(act_T)
length(annot_act_T)

````
Merge the two data sets

````
gene.in.common <- intersect( rownames(rest_T) , rownames(act_T) ) ; length(gene.in.common)
merged_raw_counts <- cbind(rest_T[gene.in.common,] , act_T[gene.in.common,] )
merged_raw_counts<-data.matrix(merged_raw_counts)
dim(merged_raw_counts)
str(merged_raw_counts)
````
Check the matrix

````
str(merged_raw_counts)
dim(merged_raw_counts)
class(merged_raw_counts)
head(merged_raw_counts[,1:3])
````
Quality controls

````
pdf("Quality_Control.pdf",onefile = F) 
par(mfrow=c(1,3))
boxplot(colSums(merged_raw_counts), main = "Sequencing depth",frame.plot=F,col="orange")
boxplot(colSums(merged_raw_counts > 0), main = "Detected genes",frame.plot=F,col="cyan")
smoothScatter(colSums(merged_raw_counts), colSums(merged_raw_counts > 0), main = "Saturation plot",frame.plot=FALSE,ylab = "Detected genes",xlab = "Sequencing depth",cex = 2,postPlotHook=NULL)
dev.off()

````

Filter for  low quality cells

````
merged_raw_counts.filtered <- merged_raw_counts[ , colSums(merged_raw_counts) > 1000 & colSums(merged_raw_counts) < 100000 ]
````
Check matrix

````
dim(merged_raw_counts.filtered)
str(merged_raw_counts.filtered)
class(merged_raw_counts.filtered)
head(merged_raw_counts.filtered[,1:3])
````

Filter for not expressed genes
````
merged_raw_counts.filtered <- merged_raw_counts.filtered[ rowSums(merged_raw_counts.filtered,na.rm =T )> 0 , ]
````
Check matrix

````
dim(merged_raw_counts.filtered)
str(merged_raw_counts.filtered)
class(merged_raw_counts.filtered)
head(merged_raw_counts.filtered[,1:3])
````

# Perepare data for clustering analysis (SCANPY pipeline in python)
````
merged_raw_counts.filtered_2<-merged_raw_counts.filtered
rownames(merged_raw_counts.filtered_2)<-substr(rownames(merged_raw_counts.filtered_2),1,15)
````
# Normalize to cpm
````
merged.cpm <- log2(t(t(merged_raw_counts.filtered_2)/(colSums(merged_raw_counts.filtered_2)/1e6)) + 1)
````

Check matrix

````
dim(merged_raw_counts.filtered_2)
str(merged_raw_counts.filtered_2)
class(merged_raw_counts.filtered_2)
head(merged_raw_counts.filtered_2[,1:3])
````
# Generate annotation files
````
annotation_D1_Lung<-substr(colnames(merged_raw_counts.filtered_2),1,5)
````
# Save the files for downstream analyses
````
write.table(merged_raw_counts.filtered_2,"~/merged_raw_counts.filtered.donor1.txt",sep="\t")
write.table(merged.cpm,"~/Normalized_merged_counts.filtered.donor1.txt",sep="\t")
write.table(as.vector(annotation_D1_Lung),"~/annotation_D1_Lung.txt",sep="\t",row.names = F,quote=F,col.names = F)
````

# Clustering analysis
There are two options: i) the first option is to perform clustering using gene expression data; ii) the second option is to apply metaVIPER for clustering analysis (using GTEX networks for normal cells and TCGA for cancer cells)

# Clustering based on GTEX networks

Set directory for GETX  networks

````
path_dir<-'~/ac_lab_scratch/CZI/aracne-repository/GTEx/'
file_names=as.list(dir(path = path_dir, pattern="*.rda"))

regList <- list()
for(i in 1:length(file_names)) 
{
  new.reg <- get(load(paste0(path_dir,file_names[[i]]))) 
  new.name <-paste0("reg_",i)
  regList[[new.name]] <- new.reg
}
````

Load the the gene expression file (normalized)

````
library(data.table)
library(atools)
exp_mat<-fread("~/Normalized_merged_counts.filtered.donor1.txt",colClasses = "numeric")
exp_mat<-as.data.frame(exp_mat)
head(exp_mat[,1:4])
ensemble_id<-exp_mat[,1]
exp_mat<-exp_mat[,-1]
rownames(exp_mat)<-ensemble_id
head(exp_mat[,1:4])

# Convert from Ensemble to entrezID, because all the GTEX networks have been generated in entrez IDs
exp_entrez<-Ensemble2entrez(exp_mat)

head(exp_entrez[,1:3])

#Compute the double rank singnature
rank_exp_entrez <- apply(exp_entrez, 2, rank)
median <- apply(rank_exp_entrez, 1, median)
mad <- apply(rank_exp_entrez, 1, mad)
signature_entrez<- (rank_exp_entrez - median)/mad

````
Apply metaVIPER

````
mat_pa<-metaVIPER(eset = signature_entrez,regulon = regList,weight = "max", method = "none") #ATTENTION: "method" must be "none" when the signature is precomputed. Otherwise metaVIPER will calculate the signature using the method "scale"

#convert to symbols
rownames(mat_pa)<-entrez2gene(rownames(mat_pa))
write.table(mat_pa,"~/PA_GTEX_D1_Lung.txt",sep="\t")
````

Clustering based on protein activity

````
library(reticulate)
command<- "python3.6"
path2script='~/ProteinActivityClusteringSCANPY.py'
args = c('PA_GTEX_D1_Lung.txt', 'annotation_D1_Lung.txt','out_Clusters_PA_GTEX_D1_LUNG.txt')
allArgs = c(path2script,args )
system2(command, allArgs)

````
## Meta cells inference based on viperSimilarity based clusters



















