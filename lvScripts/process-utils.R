#' Filters out genes with no expression and low quality cells.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param minCount Minimum number of reads in a cell. Default of 1000.
#' @param maxCount Maximum number of reads in a cell. default of 100000.
#' @return Quality controlled matrix.
QCTransform <- function(raw.mat, minCount = 1000, maxCount = 100000) {
  filt.mat <- raw.mat[, colSums(raw.mat) > minCount & colSums(raw.mat) < maxCount]
  filt.mat <- filt.mat[ rowSums(raw.mat) != 0 ,]
  rem.genes <- nrow(raw.mat) - nrow(filt.mat); rem.cells <- ncol(raw.mat) - ncol(filt.mat)
  print(rem.genes)
  print(paste('Removed ', rem.genes, ' genes and ', rem.cells, ' cells.', sep =''))
  return(filt.mat)
}

#' Generates basic quality control plots from raw gene expression data.
#' 
#' @param raw.mat Matrix of raw gene expression data (genes X samples).
#' @param plot.path Path for plots to be saved.
#' @param plot.prefix Prefix for file names of plots. If not given, no prefix will be added.
QCPlots <- function(raw.mat, plot.path, plot.prefix) {
  if (missing(plot.prefix)) {
    plot.name <- paste(plot.path, 'QC-plots.pdf', sep = '')
  } else {
    plot.name <- paste(plot.path, plot.prefix, '_QC-plots.pdf', sep = '')
  }
  ## generate plots
  pdf(plot.name, onefile = TRUE, width = 8, height = 4)
  par(mfrow=c(1,3))
  boxplot(colSums(raw.mat), main = "Sequencing depth", frame.plot=F, col="orange")
  boxplot(colSums(raw.mat > 0), main = "Detected genes", frame.plot=F, col="cyan")
  smoothScatter(colSums(raw.mat), colSums(raw.mat > 0), main = "Saturation plot",
                frame.plot=FALSE, ylab = "Detected genes", xlab = "Sequencing depth", 
                cex = 2, postPlotHook=NULL)
  dev.off()
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
