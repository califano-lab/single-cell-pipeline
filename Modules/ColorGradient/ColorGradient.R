expColor <- function(gene, exp, col = c("Grey", "Pink", "Red")){
    ind <- match(gene, rownames(exp))
    if (!is.na(ind)){
        c <- exp[ind, ]
        if (sum(c) > 0){
            c <- floor(((c-min(c))/(max(c)-min(c)))*100)+1
            c <- colorRampPalette(col)(102)[c]
        }else c <- rep(col[1], ncol(exp))
    }else c <- rep(col[1], ncol(exp))
    return(c)
}
#expression color gradient

vpColor <- function(reg, vp, col = c("Blue", "Grey", "Red"), range = 5){
    ind <- match(reg, rownames(vp))
    if (!is.na(ind)){
        c <- vp[ind, ]
        c[abs(c) > range] <- range * sign(c)[abs(c) > range]
        c <- sign(c)*(abs(c)/range)*100
        c <- sapply(c, function(x, col){
            if (x > 0) colorRampPalette(col[c(2, 3)])(102)[floor(x)+1]
            else colorRampPalette(col[c(2, 1)])(102)[floor(abs(x))+1]
        }, col=col)
    }else c <- rep(col[2], ncol(vp))
    return(c)
}
#activity color gradient

nesColor <- function(set, nes, col = c("Blue", "Grey", "Red"), range = 5){
    ind <- match(set, rownames(nes))
    if (!is.na(ind)){
        c <- nes[ind, ]
        c[abs(c) > range] <- range * sign(c)[abs(c) > range]
        c <- sign(c)*(abs(c)/range)*100
        c <- sapply(c, function(x, col){
            if (x > 0) colorRampPalette(col[c(2, 3)])(102)[floor(x)+1]
            else colorRampPalette(col[c(2, 1)])(102)[floor(abs(x))+1]
        }, col=col)
    }else c <- rep(col[2], ncol(nes))
    return(c)
}
#nES color gradient

geneColor <- function(x){
    col <- colSums(x > 0)
    col <- floor(((col-min(col))/(max(col)-min(col)))*1000)+1
    col <- colorRampPalette(c("grey", "red"))(1000)[col]
    return(col)
}
#gene number color gradient

countColor <- function(count){
    col <- log2(colSums(count))
    col <- floor(((col-min(col))/(max(col)-min(col)))*1000)+1
    col <- colorRampPalette(c("grey", "red"))(1000)[col]
    return(col)
}
#count number color gradient
