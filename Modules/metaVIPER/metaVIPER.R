metaVIPER <- function(eset, regulon, weight = c("mean", "max", "nes", "user"), ts = NA,
                      dnull = NULL, pleiotropy = FALSE, nes = TRUE,
                      method = c("scale", "rank", "mad", "ttest", "none"), bootstraps = 0,
                      minsize = 25, adaptive.size = FALSE, eset.filter = TRUE,
                      pleiotropyArgs = list(regulators = 0.05, shadow = 0.05, targets = 10,
                                            penalty = 20, method = "adaptive"),
                      cores = 1, verbose = TRUE){
  if (is(eset, "viperSignature")){
    ns <- list(ncol(eset$signature), colnames(eset$signature))
  }else if (is(eset, "matrix")){
    ns <- list(ncol(eset), colnames(eset))
  }#object check of eset
  regulators <- unique(unlist(lapply(regulon, function(x) names(x))))
  #regulators used in the analysis
  
  message("\nExtracting regulatory information")
  vp <- lapply(regulon, function(x, eset, dnull, pleiotropy, nes, method, bootstraps,
                                 minsize, adaptive.size, eset.filter, pleiotropyArgs, cores, verbose, ns, regulators){
    vp <- matrix(0, length(regulators), ns[[1]], dimnames = list(regulators, ns[[2]]))
    vpmat <- viper(eset = eset, regulon = x, dnull = dnull, pleiotropy = pleiotropy, nes = nes,
                   method = method, bootstraps = bootstraps, minsize = minsize, adaptive.size = adaptive.size,
                   eset.filter = eset.filter, pleiotropyArgs = pleiotropyArgs, cores = cores, verbose = verbose)
    vp[rownames(vpmat), ] <- vpmat
    vp
  }, eset = eset, dnull = dnull, pleiotropy = pleiotropy, nes = nes, method = method, bootstraps = bootstraps,
  minsize = minsize, adaptive.size = adaptive.size, eset.filter = eset.filter,
  pleiotropyArgs = pleiotropyArgs, cores = cores, verbose = verbose, ns = ns, regulators = regulators)
  #compute viper activity using each interactome
  
  message("\nIntegrating regulatory information\n")
  vp <- array(unlist(vp), dim = c(length(regulators), ns[[1]], length(regulon)),
              dimnames = list(regulators, ns[[2]], names(regulon)))
  if (weight == "mean"){
    wm <- apply(vp, c(1, 2), mean)
  }else if (weight == "max"){
    wm <- apply(vp, c(1, 2), function(x) max(abs(x)))
  }else if (weight == "nes"){
    w <- apply(abs(vp), c(1, 2), sum)
    ws <- apply(vp*abs(vp), c(1, 2), sum)
    wm <- ws/w
  }else if (weight == "user"){
    if (is.na(ts)){
      message("Error, no weight specified")
    }else{
      wm <- rowSums(sapply(1:length(regulon), function(i) vp[, , i]*ts[i]))
      wm <- matrix(wm, nrow = length(regulators), ncol = ns[[1]],
                   dimnames = list(regulators, ns[[2]]))
    }
  }
  #integrating viper activity produced by individual ARACNe network
  return(wm)
}
