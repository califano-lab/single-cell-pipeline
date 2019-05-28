#' Filters out genes with no expression and low quality cells.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of reads in a cell. Default of 1000.
#' @param maxCount Maximum number of reads in a cell. Default of 100000.
#' @param minGeneReads Minimum number of reads for a gene to be kept. Default of 1 (any gene with no reads will be removed).
#' @return Quality controlled matrix.
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000, minGeneReads = 1) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(raw.mat) >= minGeneReads ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(paste('Removed ', rem.genes, ' genes and ', rem.cells, ' cells.', sep =''))
  return(filt.mat)
}

#' Generates basic quality control plots from raw gene expression data.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param plot.path Optional argumetn of save path for plot.'
QCPlots <- function(raw.mat, plot.path) {
  if (!missing(plot.path)) {
    pdf(plot.path)
  }
  ## generate plots
  nf <- layout(matrix(c(1,2,3), nrow = 1), widths = c(3, 3, 3), heights = c(5, 5, 5), TRUE)
  boxplot(colSums(raw.mat), main = "Sequencing depth", frame.plot=F, col="orange")
  boxplot(colSums(raw.mat > 0), main = "Detected genes", frame.plot=F, col="cyan")
  smoothScatter(colSums(raw.mat), colSums(raw.mat > 0), main = "Saturation plot",
                frame.plot=FALSE, ylab = "Detected genes", xlab = "Sequencing depth", 
                cex = 2, postPlotHook=NULL)
  if (!missing(plot.path)) {
    dev.off()
  }
}

#' Performas a CPM normalization on the given data. 
#' 
#' @param dat.mat Matrix of gene expression data (genes X samples).
#' @param l2 Optional log2 normalization switch. Default of False.
#' @return Returns CPM normalized matrix
CPMTransform <- function(dat.mat, l2 = FALSE) {
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}

#' Performs a rank transformation on a given matrix.
#' 
#' @param dat.mat Matrix of data, usually gene expression (genes X samples).
#' @return Rank transformed matrix.
RankTransform <- function(dat.mat) {
  rank.mat <- apply(dat.mat, 2, rank)
  median <- apply(rank.mat, 1, median)
  mad <- apply(rank.mat, 1, mad)
  rank.mat <- (rank.mat - median) / mad
  return(rank.mat)
}

#' Transforms gene names from Ensembl to hgnc.
#' 
#' @param dat.mat Matrix of data with ENSEMBL names (genes X samples).
#' @return Data with HGNC names. Some data will likely be lost in the conversion.
Ensemble2GeneName<-function(dat.mat) {
  # packages
  require(biomaRt)
  # get the mart
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  name.map <- getBM(attributes=c('hgnc_symbol','hgnc_id','ensembl_gene_id'), 
                       filters = 'ensembl_gene_id', values = (rownames(dat.mat)), mart = ensembl)
  rownames(name.map) <- make.unique(name.map$ensembl_gene_id)
  # create a converted dataset
  convert.dat <- merge(name.map, dat.mat, by = c('row.names'))
  rownames(convert.dat) <- make.unique(convert.dat$hgnc_symbol)
  convert.dat <- convert.dat[,-c(1:4)]
  return(as.matrix(convert.dat))
}

#' Saves a matrix in a format for input to ARACNe
#'
#' @param dat.mat Matrix of data (genes X samples).
#' @param out.file Output file where matrix will be saved.
#' @param subset Switch for subsetting the matrix to 500 samples. Default TRUE.
ARACNeTable <- function(dat.mat, out.file, subset = TRUE) {
  dat.mat <- dat.mat[!duplicated(rownames(dat.mat)), ]
  if (subset) {
    dat.mat <- dat.mat[, sample(colnames(dat.mat), 500) ]
  }
  sample.names <- colnames(dat.mat)
  gene.ids <- rownames(dat.mat)
  m <- dat.mat
  mm <- rbind( c("gene", sample.names), cbind(gene.ids, m))
  write.table( x = mm , file = out.file , 
               sep="\t", quote = F , row.names = F , col.names = F )
}

#' Builds a UMAP transformation based on the shared, cell-by-cell master regulators. Utilizes custom UMAP parameters ('pearson', 25 neighbors).
#'
#' @param dat.mat Matrix of protein activity (protein X samples).
#' @param num.mrs Number of top master regulators to take from each cell. Default of 50.
#' @return Returns a UMAP based on the cell-by-cell master regulators.
cbcMR_UMAP <- function(dat.mat, num.mrs = 50) {
  require(umap)
  # set UMAP parameters
  umap_custom <- umap.defaults
  umap_custom$n_neighbors <- 25
  umap_custom$metric <- 'pearson'
  # identify MRs
  cbc.mrs <- apply(dat.mat, 2, function(x) { names(sort(x, decreasing = TRUE))[1:num.mrs] })
  cbc.mrs <- unique(unlist(as.list(cbc.mrs)))
  # generate and return UMAP
  cbc.umap <- umap(t(dat.mat[ match(cbc.mrs, rownames(dat.mat)) , ]), config = umap_custom)
  return(cbc.umap)
}

#' Unwraps a nested MR list: previous functions return cluster specific master regulators as a list of lists. This funciton will unwrap that object into one, unique list.
#'
#' @param MRs List of lists, with MR names as sub-list names and MR activity as sub-list entries.
#' @return Returns a de-duplicated list of MRs.
MR_UnWrap <- function(MRs) {
  return( unique(unlist(lapply(MRs, names), use.names = FALSE)) )
}