

# Paths for the HPC for datasets
* Datasets: 
  * ```/ifs/scratch/c2b2/ac_lab/CZI/peter_data```
  
* ARACNe repository: 
  * ```/ifs/scratch/c2b2/ac_lab/CZI/aracne-repository```

# Example DONOR1 (Lung)

````
rm(list=ls(all=T))
library(cluster)
library(ggfortify)

````
### Load data sets

Load "resting T cells" matrix

````
rest_T<-read.table("~/PP001swap.filtered.matrix.txt.bz2",sep="\t")
ensemble_id<-rest_T[,1]
dim(rest_T)
str(rest_T)
head(rest_T[,1:4])
rest_T<-rest_T[,-c(1,2)]
rownames(rest_T)<-ensemble_id
colnames(rest_T)<-rep("PP001",length(colnames(rest_T)))
head(rest_T[,1:4])
````
Load  "activated T cells" matrix

````
act_T<-read.table("~/ac_lab_scratch/CZI/peter_data/PP002swap.filtered.matrix.txt.bz2",sep="\t")
dim(act_T)
str(act_T)
ensemble_id<-act_T[,1]
act_T<-act_T[,-c(1,2)]
rownames(act_T)<-ensemble_id
head(act_T[,1:4])
colnames(act_T)<-rep("PP002",length(colnames(act_T)))
head(act_T[,1:4])
````
Merge the two data sets

````
gene.in.common <- intersect( rownames(rest_T) , rownames(act_T) ) ; length(gene.in.common)
merged_raw_counts <- cbind( rest_T[gene.in.common,] , act_T[gene.in.common,] )
dim(merged_raw_counts)
str(merged_raw_counts)
````


Filter for  low quality cells
````
merged_raw_counts.filtered <- merged_raw_counts[ , colSums(merged_raw_counts) > 1000 & colSums(merged_raw_counts) < 100000 ]
dim(merged_raw_counts.filtered)
str(merged_raw_counts.filtered)
````

Filter for not expressed genes
````
merged_raw_counts.filtered <- merged_raw_counts.filtered[ rowSums(merged_raw_counts.filtered) != 0 , ]
str(merged_raw_counts.filtered)
dim(merged_raw_counts.filtered)
````
# Normalize to cpm
merged.cpm <- log2(t(t(merged_raw_counts.filtered)/(colSums(merged_raw_counts.filtered)/1e6)) + 1)

# Clustering analysis

Select top 500 variant genes
````
var_filt<-apply(merged.cpm,1, var)
var_filt_genes <- names(sort(var_filt,decreasing = T)[1:500])
cpm_var<-merged.cpm[var_filt_genes,]
````

Define optimal number of clusters based on silhouette score
````
sil_score_ge<-NULL

#create dissimilarity matrix
diss_m<-as.matrix(1-cor(cpm_var))

for (i in 2:6)
{
  sil_score_ge<-c(sil_score_ge,pam(diss_m,i)$silinfo$avg.width)
}
# Plot the the silhouette score values and identify the optimal number of cluster
plot(c(2:6),sil_score_ge,type = "b")  
# The optimal number of cluster is 2, however these two clusters are weak
clusters_k2<-pam(diss_m,2)
````
# Principal component analysis
````
pca<-prcomp(t(cpm_var))

#Plot the PCA 
autoplot(pca)
````
# Select cells from each cluster and prepare data for network analysis
Ideally, we want to generate an ARACNe network  for each cluster by randolmly selecting "n" cells (e.g. n=1000) from each cluster. 
````
names_cells<-c(sample(names(clusters_k2$clustering)[grep("1",clusters_k2$clustering)],1000),sample(names(clusters_k2$clustering)[grep("2",clusters_k2$clustering)],1000))

# Alternatively, in case the clusters are not robust (not significant), we can randomly select cells from the entire matrix
# NOT RUN
#names_cells<-sample(colnames(cpm_var), 1000)

# check how many cells are activated T cells and how many are  resting T cells 

length(grep("PP001", names_cells)) # n=478
length(grep("PP002", names_cells)) #n=522
````

## Prepare expression matrix for ARACNe
````
expmat<-merged.cpm[,names_cells]
dim(expmat)
expmat<-unique(expmat)
dim(expmat)
#remove the version of the ensembleID: remove the ".N"
rownames(expmat)<-substr(rownames(expmat),1,15)
dim(unique(expmat))

#save the matrix 
save(expmat,file="~/ac_lab_scratch/CZI/Pas_results/DONOR1/LUNG/expression4ARACNe.rda")
````












