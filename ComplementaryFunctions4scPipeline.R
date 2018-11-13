gene2entrez<-function (gene, sp = "Homo")
{

  if (!exists("desc", envir = as.environment(".GlobalEnv")))
    data(desc, envir = as.environment(".GlobalEnv"))
  desc <- get("desc", env = as.environment(".GlobalEnv"))
  desc <- desc[desc[, 1] %in% which.specie(sp), ]
  res <- as.numeric(desc[match(tolower(gene), tolower(desc[,3])), 2])
  names(res) <- gene
  res
}


entrez2gene<-function (id)
{
  if (!exists("desc", envir = as.environment(".GlobalEnv")))
    data(desc, envir = as.environment(".GlobalEnv"))
  desc <- get("desc", env = as.environment(".GlobalEnv"))
  res <- desc[match(id, desc[, 2]), 3]
  names(res) <- id
  res
}


competitive.test2<-function (Pvalue, Weight, p_random = NA) 
{
  if (is.na(p_random)[1]) {
    p_random <- matrix(runif(length(Pvalue) * 1e+03), nrow = length(Pvalue))
  }
  if (nrow(p_random) != length(Pvalue)) {
    stop("Error: nrow(p_random)!=length(Pvalue)")
  }
  random <- NULL
  for (i in 1:ncol(p_random)) {
    random[i] <- selfcontained.test(pvalue = p_random[, i], 
                                    weight = Weight, p_permu = NA)[[1]]
  }
  return(list(`significance level for combining pvalues` = mean(selfcontained.test(Pvalue,Weight, p_permu = NA)[[1]] > random)))
}


stouff<-function(NES)
{
  sum(NES)/sqrt(length(NES))
}





Ensemble2GeneName<-function(dataset2Convert)
{
require(biomaRt)
#ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
 load("/Users/pl2659/Documents/GitHub/TME_PDA/input_data/ensemblMay2017.rda")
  names_dataset<-getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'),filters = 'ensembl_gene_id', values = (rownames(dataset2Convert)), mart = ensemblMay2017)
  rownames(names_dataset)<-make.unique(names_dataset$ensembl_gene_id)
  dataset2Convert_2<-merge(names_dataset,dataset2Convert,by=c("row.names"))
  dim(dataset2Convert_2)
  rownames(dataset2Convert_2)<-make.unique(dataset2Convert_2$hgnc_symbol)
  GeneName_dataset<-dataset2Convert_2[,-c(1:4)]
  head(GeneName_dataset[,1:4])
  return(GeneName_dataset)
}


Ensemble2entrez<-function(dataset2Convert)
{
  require(biomaRt)
  #ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  load("/Users/pl2659/Documents/GitHub/TME_PDA/input_data/ensemblMay2017.rda")
  names_dataset<-getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id',"entrezgene"),filters = 'ensembl_gene_id', values = (rownames(dataset2Convert)), mart = ensemblMay2017)
  rownames(names_dataset)<-make.unique(names_dataset$ensembl_gene_id)
  dataset2Convert_2<-merge(names_dataset,dataset2Convert,by=c("row.names"))
  dim(dataset2Convert_2)
  entrez_unique<-make.names(dataset2Convert_2$entrezgene,unique = T)
  rownames(dataset2Convert_2)<-substring(entrez_unique, 2)
  GeneName_dataset<-dataset2Convert_2[,-c(1:5)]
  head(GeneName_dataset[,1:5])
  return(GeneName_dataset)
}



UNIPROT2entrezID<-function(vector2Convert)
{
  require(AnnotationDbi)
  require(org.Hs.eg.db)
  cols <- c("SYMBOL", "ENTREZID")
  uniKeys<-as.character(vector2Convert)
  entrez_ids<-select(org.Hs.eg.db, keys=uniKeys, columns=cols, keytype="UNIPROT",multiVals="asNA")
  return(entrez_ids)
}




gsea_Clust2<-function(signature_gene,list_sets)
{
 require(atools)
 signature_gene<-na.omit(signature_gene)
  mat0<-matrix(nrow = length(list_sets),ncol = length(signature_gene[1,]))
  rownames(mat0)<-names(list_sets)
  colnames(mat0)<-colnames(signature_gene)
  for (i in 1:length(list_sets) )
  {
    mat0[i,]<-apply(signature_gene,2,function(x){res<-gsea(c(x),list_sets[[i]],per=100,pout = F)$nes})
  }
  return(mat0)
}


colHeat = function(gene, rpm, colFrom="Blue", colTo="Red", colGrad=50){
				   ind = rpm[match(gene, rownames(rpm)), ]
				   ind = floor(((ind-min(ind))/(max(ind)-min(ind)))*colGrad)+1
				   return(colorRampPalette(c(colFrom, colTo))(colGrad)[ind])
}		



Knn_metaCells<-function(normData, rawCounts, K)
{
  require(FNN)
  #Get the KNN using the clustering onnormalized expression data
  knn_10<-get.knn(t(normData),K,algorithm = "brute") 
  
  # Create a matrix of the same size of expmat0
  sum_knn<-matrix(0,length(normData[,1]),length(normData[1,]))
  mdata_c1<-rawCounts[,colnames(normData)]
  
  message("Check size")
  dim(mdata_c1)
  
  ## check that  the same  cells are selected
  for (i in 1:length(mdata_c1[1,]))
  {
    sum_knn[,i]<-apply(mdata_c1[,colnames(mdata_c1[,c(i,knn_10$nn.index[i,])])],1,sum)
    print(paste0("MetaCell_",i))
  }
  dim(sum_knn)
  rownames(sum_knn)<-rownames(mdata_c1)
  colnames(sum_knn)<-colnames(mdata_c1)
  
  ind <- colSums(sum_knn)>0
  ### Now let's compute tpm, signature and protein activity
  tpm_KNN <- log2(t(t(sum_knn[, ind])/(colSums(sum_knn[, ind])/1e6)) + 1)
  return(tpm_KNN)
}


getKNN_VIPER<- function(ViperMatrix,K)
  {
  #create an empty matrix
  dat_knn<-matrix(0,length(ViperMatrix[,1]),K)
  for(i in 1:length(ViperMatrix[,1]))
  {
  dat_knn[i,]<-sort(ViperMatrix[2,],decreasing = T,index.return=TRUE)$ix[1:K]
  
  return(dat_knn)
  }