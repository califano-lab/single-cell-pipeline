args <- commandArgs(trailingOnly = T)
wd <- args[1]#working directory
expfn <- args[2]#expression file
regfn <- args[3]#regulators
mht <- args[4]#correction (none, bonferroni, fdr)
alpha <- as.numeric(args[5])#p-value
#arguments

fn <- list.files(path = wd, pattern = "bootstrapNetwork", full.names = T)
res <- resmi <- edges <- NULL
for (fn1 in fn) {
    message(fn1)
    tmp <- strsplit(readLines(fn1), "\t")
    mi <- as.numeric(sapply(tmp, function(x) x[3]))
    names(mi) <- sapply(tmp, function(x) paste(x[1:2], collapse="--"))
    edges <- c(edges, length(mi))
    ed <- unique(c(names(res), names(mi)))
    pos1 <- match(ed, names(res))
    pos1na <- which(is.na(pos1))
    pos2 <- match(ed, names(mi))
    pos2na <- which(is.na(pos2))
    res <- res[pos1]
    resmi <- resmi[pos1]
    res[pos1na] <- resmi[pos1na] <- 0
    names(res)[pos1na] <- names(resmi)[pos1na] <- ed[pos1na]
    mi <- mi[pos2]
    mi[pos2na] <- 0
    res <- res+sign(mi)
    resmi <- resmi + mi
}
resmi <- resmi/res
regs <- readLines(regfn)
genes <- sapply(strsplit(readLines(expfn)[-1], "\t"), function(x) x[1])
tedges <- length(genes)*length(which(genes %in% regs))
lambda <- sum(edges)/tedges
res <- p.adjust(ppois(res, lambda, lower.tail = F), mht)
res <- sort(res[res < alpha])
tmp <- cbind(t(sapply(strsplit(names(res), "--"), function(x) x)), signif(resmi[match(names(res), names(resmi))], 5), signif(res, 5))
cat(unlist(apply(tmp, 1, paste, collapse="\t"), use.names=FALSE), sep="\n", file = file.path(wd, "finalNetwork_4col.tsv"))
message("Finished!")
#consolidate bootstrap networks
