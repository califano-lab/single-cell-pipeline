

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

i) Gene expression based clustering: APPLY LOUVEIN CLUSTERING ALGORITHM (from SCANPY pipeline)
````
library(reticulate)

command<- "python3.6"
path2script='~/Cluster_Exp.py'
args = c('merged_raw_counts.filtered.donor1.txt', 'annotation_D1_Lung.txt','out_Clusters_D1_LUNG.txt')
allArgs = c(path2script,args )
system2(command, allArgs)
````

Otherwise go to the terminal and type:
````
python3.6 Cluster_Exp.py merged_raw_counts.filtered.donor1.txt  annotation_D1_Lung.txt out_Clusters_D1_LUNG.txt
````

Prepare data for ARACNe
 ````
clusters<-read.table("/Users/pl2659/Documents/SCANPY/out_Clusters_D1_LUNG.txt",sep="\t",header = T)
head(clusters)
message("Check the number of clusters")
max(clusters$louvain)+1 # python start from zero, the if the max is 4 the number of clusters are 5 

cluster_cells<-NULL
for( i in 0:max(clusters$louvain))
{
  cluster_cells<-c(cluster_cells,sum((clusters$louvain==i)=="TRUE"))
  message(paste("Cluster",i,'= '),sum((clusters$louvain==i)=="TRUE"))
}

message("Select all the cells with more than 300 clusters")
names(cluster_cells)<-c(0:max(clusters$louvain))
barplot(sort(cluster_cells),horiz = T)
ind_clusters<-names(which(cluster_cells>300))

cluster_0_cells<-clusters$X[grep("0",clusters$louvain)]
cluster_1_cells<-clusters$X[grep("1",clusters$louvain)]
cluster_2_cells<-clusters$X[grep("2",clusters$louvain)]

expmat0<-merged.cpm[,cluster_0_cells]
expmat1<-merged.cpm[,cluster_1_cells]
expmat2<-merged.cpm[,cluster_2_cells]

dim(expmat0)
dim(expmat1)
dim(expmat2)

````
Save the matrices for ARACNe

````
saveRDS(expmat0,file="~/d1-lung_c0_expression4ARACNe.rds")
saveRDS(expmat1,file="~/d1-lung_c1_expression4ARACNe.rds")
saveRDS(expmat2,file="~/d1-lung_c2_expression4ARACNe.rds")
````








Meta-cells inference in each cluster
````
library(data.table)
library(reticulate)
library(atools)
source('~ComplementaryFunctions4scPipeline.r')
````
load the  the raw count UMI 
````
mdata<-fread("~/merged_raw_counts.filtered.donor1.txt",colClasses = "numeric")

mdata<-as.data.frame(mdata)
ensemble_id<-mdata[,1]
mdata<-mdata[,-1]
rownames(mdata)<-ensemble_id
````

# Build meta-cells for each cluster using the function "Knn_metaCells"

````
meta_expmat0<-Knn_metaCells(expmat0,mdata,10)
meta_expmat1<-Knn_metaCells(expmat1,mdata,10)
meta_expmat2<-Knn_metaCells(expmat2,mdata,10)
head(meta_expmat0[,1:4])
head(tpm_KNN[,1:4])

## check size
dim(meta_expmat0)
dim(meta_expmat1)
dim(meta_expmat2)
#Save the matrix 
save(meta_expmat0,file="/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/MC_d1-lung_c0_expression4ARACNe.rda")
save(meta_expmat1,file="/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/MC_d1-lung_c1_expression4ARACNe.rda")
save(meta_expmat2,file="/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/MC_d1-lung_c2_expression4ARACNe.rda")

````
Then, you need to randomly  select >200 cells  for each cluster and proceed with ARACNe to build a network for each cluster.
The same procedure is applied if clustering is performed using metaVIPER.

# Generate ARACNe networks for each cluster (with more than 300 cells)

Please visit  the ARACNe repository before to run ARACNe
````
sh ARACNe_p.sh
````
# Virtual inference of protein activity  analysis by metaVIPER
````
library(viper)
source("R/ComplementaryFunctions.r")
load("R/desc.rda")

eset<-ExpressionSet(expmat)
# Example  of hoe to generate a regulon object
net_path<-file.path("~/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/LUNG_DONOR1-tf-network.tsv")
example_reg<-aracne2regulon(net_path,eset,format = "3col")
#save(example_reg,file="~/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/regul_LUNG_DONOR.rda")

After all the regulons are computed we can load them and apply metaVIPER

#Load the regulons
net1<-readRDS("/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/d1-lung_cluster0_mergedNet.rds")
net2<-readRDS("/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/d1-lung_cluster1_mergedNet.rds")
net3<-readRDS("/Volumes/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/d1-lung_cluster2_mergedNet.rds")

## Load expression matrix (normalized)for each cluster

exp_cluster1<-readRDS(~/d1-lung_c0_expression4ARACNe.rds)
# Infer  protein activity 
pa_D1_lung<-metaVIPER(merged.exp_cluster1, regulon =as.list(net1,net2,net3),method = "scale") #*

## write the Protein activity matrix and do clustering analysis
write.table(pa_D1_lung,"~/PA__D1_lung_metaVIPER.txt",sep="\t")
````
*You can also use the double-rank transformation to compute signature but, in this case, remember to set method="none".

This is the code to run the double rank signature
````
# DON'T RUN 
#rank_exp <- apply(exp_cluster1, 2, rank)
#median_exp <- apply(rank_exp, 1, median)
#mad_exp <- apply(rank_exp, 1, mad)
#signature_exp <- (rank_exp - median_exp)/mad_exp
````

# Clustering analysis based on protein activity (python script: "ProteinActivityClusteringSCANPY.py")

````
library(reticulate)
command<- "python3.6"
path2script='~/ProteinActivityClusteringSCANPY.py'
args = c('PA__D1_lung_metaVIPER.txt', 'annotation_D1_Lung.txt','out_Clusters_PA_D1_LUNG.txt')
allArgs = c(path2script,args )
system2(command, allArgs)
````
#  Select cells associated to each cluster and indentify the MRs using the Stouffer intergration.


















