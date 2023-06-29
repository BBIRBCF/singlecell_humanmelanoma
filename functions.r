prepareGS <- function(gsets, GENES, minGenes=10, maxGenes=2000){
    GENESETS <- lapply(1:length(gsets), function(S){
        x <- readLines(gsets[S])
        AUX <- sapply(x, function(y){
            aux <- strsplit(y,"\t")[[1]]
               geneset <- as.character(aux[c(1,3:length(aux))])
            return(geneset)
        })
        names(AUX) <- sapply(AUX, function(x) x[1])
        AUX
    })
    names(GENESETS) <- basename(gsets)

    GENESETS2 <-list()
    for(i in 1:length(gsets)){
        GENESETS2[[i]] <- list()
        for(j in 1:length(GENESETS[[i]])){
            x <- GENESETS[[i]][[j]]; x <- c(x[1], x[x %in% GENES])
            if(length(x)<=maxGenes-1&length(x)>=minGenes+1)
                GENESETS2[[i]][[ x[1] ]] <- x[2:length(x)]
        }
    }
    names(GENESETS2) <- basename(gsets)
    GENESETS2
}


runHyper <- function(GENESETS2, gsets, forD, res){
    overlaps <- list()
    for(o in names(forD)){
        int.geneset <- forD[[o]]
        overlaps[[o]] <- lapply(1:length(GENESETS2), function(S){
            RET.LIST <- t(sapply(1:length(GENESETS2[[S]]),function(k){
                geneset <- GENESETS2[[S]][[k]]
                ngs <- length(geneset)
            nsig <- length(int.geneset)
                nsig.gs <- length(intersect(geneset,int.geneset))
                npop <- length(GENES)
                pv <- phyper(nsig.gs-1, ngs, npop-ngs, nsig, lower.tail=F)
                c(names(GENESETS2[[S]])[k], ngs, nsig, nsig.gs, npop, pv, paste(intersect(geneset,int.geneset),collapse=" ") )
            }))
            colnames(RET.LIST) <- c("path","n.geneset.path", "n.geneset.int","n.intersect", "n.population", "p.val", "common.genes")
            RET.LIST <- as.data.frame(RET.LIST)
            RET.LIST$p.val <- as.numeric(as.character(RET.LIST$p.val))
            RET.LIST$p.val.adj <- p.adjust(RET.LIST$p.val, method="BH")
            RET.LIST <- RET.LIST[order(RET.LIST$p.val), ]
            RET.LIST <- RET.LIST[,c("path", "n.geneset.path", "n.geneset.int","n.intersect", "n.population", "p.val", "p.val.adj", "common.genes")]
            write.table(RET.LIST, file= paste0(res, o,".",basename(gsets)[S],".xls"),row.names=FALSE, quote=FALSE,sep="\t")
            RET.LIST
        })
        names(overlaps[[o]]) <- basename(gsets)
    }
    overlaps
}

runHyperhtml  <-
function(GENESETS2, gsets, forD, res){
    overlaps <- list()
    for(o in names(forD)){
        int.geneset <- forD[[o]]
        overlaps[[o]] <- lapply(1:length(GENESETS2), function(S){
            RET.LIST <- t(sapply(1:length(GENESETS2[[S]]),function(k){
                geneset <- GENESETS2[[S]][[k]]
                ngs <- length(geneset)
            nsig <- length(int.geneset)
                nsig.gs <- length(intersect(geneset,int.geneset))
                npop <- length(GENES)
                pv <- phyper(nsig.gs-1, ngs, npop-ngs, nsig, lower.tail=F)
                c(names(GENESETS2[[S]])[k], ngs, nsig, nsig.gs, npop, pv, paste(intersect(geneset,int.geneset),collapse=" ") )
            }))
            colnames(RET.LIST) <- c("path","n.geneset.path", "n.geneset.int","n.intersect", "n.population", "p.val", "common.genes")
            RET.LIST <- as.data.frame(RET.LIST)
            RET.LIST$p.val <- as.numeric(as.character(RET.LIST$p.val))
            RET.LIST$p.val.adj <- p.adjust(RET.LIST$p.val, method="BH")
            RET.LIST <- RET.LIST[order(RET.LIST$p.val), ]
            RET.LIST <- RET.LIST[,c("path", "n.geneset.path", "n.geneset.int","n.intersect", "n.population", "p.val", "p.val.adj", "common.genes")]
            p <- openPage(  paste0(res, o,".",basename(gsets)[S],".html"))
            hwrite(RET.LIST, p)
            closePage(p)
            RET.LIST
        })
        names(overlaps[[o]]) <- basename(gsets)
    }
    overlaps
}


buildResHTML <- function(res, forD){
    resdirs <- dir(res)
    resdirs <- resdirs[regexpr(".xls",resdirs)>0]
	resdirs <- resdirs[order(factor(sapply(strsplit(resdirs, ".", fixed=TRUE),"[",1), levels = names(forD)))]
    resdirs <- lapply(0:(length(forD)-1), function(i) resdirs[((i*6)+1):((i+1)*6)])
    resdirs <- do.call(cbind, resdirs)
    rownames(resdirs) <- c("Broad Hallmarks", "GOBP", "GOCC", "GOMF", "GOSLIM", "KEGG")
    colnames(resdirs) <- sub("tive.*", "tive", sub(res, "", resdirs[1,]))
    xout <- matrix("link", ncol=ncol(resdirs), nrow=nrow(resdirs))
    for(i in 1:ncol(resdirs)) xout[,i] <- rownames(resdirs)
    colnames(xout) <- colnames(resdirs); rownames(xout) <- rownames(resdirs)
    xout <- data.frame(xout)
    links <- lapply(1:ncol(resdirs), function(i) resdirs[,i])
    list(xout=xout, links=links)
}
