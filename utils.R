#################################################################################################################
#################################################################################################################
### utils.R: useful functions
### Author: Toni Berenguer
#################################################################################################################
#################################################################################################################

plotAffyRNAdeg.mod <- function (rna.deg.obj, transform = "shift.scale", cols = NULL, ...) 
{
    if (!is.element(transform, c("shift.scale", "shift.only", 
        "neither"))) 
        stop("Tranform must be 'shift.scale','shift.only', or 'neither'")
    mns <- rna.deg.obj$means.by.number
    if (is.null(cols)) 
        cols = rep(4, dim(mns)[1])
    ylab = "Mean Intensity"
    if (transform == "shift.scale") {
        sds <- rna.deg.obj$ses
        mn <- mns[, 1]
        mns <- sweep(mns, 1, mn)
        mns <- mns/(sds)
        mns <- sweep(mns, 1, 1:(dim(mns)[1]), "+")
        ylab <- paste(ylab, ": shifted and scaled")
    }
    else if (transform == "shift.only") {
        mn <- mns[, 1]
        mns <- sweep(mns, 1, mn)
        mns <- sweep(mns, 1, 1:(dim(mns)[1]), "+")
        ylab <- paste(ylab, ": shifted")
    }
    plot(-2, -1, pch = "", xlim = range(-1, (dim(mns)[2])), ylim = range(min(as.vector(mns)) - 
        1, max(as.vector(mns)) + 1), xlab = "5' <-----> 3'\n Probe Number ", 
        ylab = ylab, axes = FALSE, main = "RNA degradation plot", 
        ...)
    axis(1)
    axis(2)
    for (i in 1:dim(mns)[1]) lines(0:((dim(mns)[2] - 1)), mns[i, 
        ], col = cols[i])
	invisible(mns)
}

###################################################################################################################
    # write.html functio from Evarist modified to allow center column option
    
write.html.mod <- function (x, links, tiny.pic, tiny.pic.size = 100, title = "",
    file, digits = 3, col.align='center', cellpadding=10)
{
    stopifnot(class(x) == "data.frame")
    if (missing(links))
        links <- vector("list", ncol(x))
    if (missing(tiny.pic))
        tiny.pic <- vector("list", ncol(x))
    stopifnot(class(links) == "list")
    stopifnot(class(tiny.pic) == "list")
    stopifnot(length(links) == ncol(x))
    stopifnot(length(tiny.pic) == ncol(x))
    stopifnot(!missing(file))
    column.class <- unlist(lapply(x, class))
    for (j in 1:ncol(x)) {
        if (column.class[j] == "factor")
            x[, j] <- as.character(x[, j])
        if (column.class[j] == "numeric")
            x[, j] <- round(x[, j], digits = digits)
    }
    cat("<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n",
        sep = "", file = file)
    cat("<html>\n", file = file, append = T)
    cat("<body>\n", file = file, append = T)
    cat(paste("<CAPTION ALIGN=\"top\"><center><B>", title, "</B></center></CAPTION><BR>\n"),
        sep = "", file = file, append = T)
    cat(paste("<TABLE border=1 cellpadding=", cellpadding,">\n", sep=''), file = file, append = T)
    cat("<TR>\n", file = file, append = T)
    for (j in 1:ncol(x)) {
        cat("<TH>", file = file, append = T)
        cat(colnames(x)[j], file = file, append = T)
        cat("</TH>\n", file = file, append = T)
    }
    cat("</TR>\n", file = file, append = T)
    for (i in 1:nrow(x)) {
        cat("<TR>\n", file = file, append = T)
        for (j in 1:ncol(x)) {
            cat(paste("<TD align=", col.align, ">", sep=''), file = file, append = T)
            if (is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(x[i, j], file = file, append = T)
            }
            else if (is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                  tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                  "\" width=\"", tiny.pic.size, "\" /></A>",
                  sep = ""), file = file, append = T)
            }
            else if (!is.null(links[[j]]) & is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\">",
                  x[i, j], "</A>", sep = ""), file = file, append = T)
            }
            else if (!is.null(links[[j]]) & !is.null(tiny.pic[[j]])) {
                cat(paste("<A HREF=\"", links[[j]][[i]], "\"><img src=\"",
                  tiny.pic[[j]][[i]], "\" height=\"", tiny.pic.size,
                  "\" width=\"", tiny.pic.size, "\" /></A>",
                  sep = ""), file = file, append = T)
            }
            cat("</TD>\n", file = file, append = T)
        }
        cat("</TR>\n", file = file, append = T)
    }
    cat("</TABLE>\n", file = file, append = T)
    cat("</body>\n", file = file, append = T)
    cat("</html>\n", file = file, append = T)
    sortDragHtmlTable(filename = file)
}

#################################################################################################################
	# QQ-plot for p-values


qqpval<-function(p, pch=16, col=4, ...)
{
   p<-p[!is.na(p)]
   n<-length(p)
   pexp<-(1:n)/(n+1)
   plot(-log(pexp,10), -log(sort(p),10), xlab="-log(expected P value)",
        ylab="-log(observed P value)", pch=pch, col=col, ...)
   abline(0,1,col=2)
}

#################################################################################################################

make.contrasts <- function(d, form, lsel, signs, type='balanced', ngroups=NULL)
{
   
    nv <- strsplit(paste(form)[2], split=' \\+ ')[[1]];
    nv <- nv[regexpr("\\:", nv) < 0];
    nv <- nv[nv!='-1']
    for (n in nv) if (is.character(d[, n])) d[, n] <- factor(d[, n])

    nvf <- nv[sapply(nv, function(v, d) is.factor(d[, v]), d)];
    nvn <- nv[!nv%in%nvf];

    daf <- expand.grid(sapply(nvf, function(o) levels(d[, o]), simplify=F));
    dan <- t(matrix(apply(d[, nvn, drop=F], 2, mean), ncol=nrow(daf)));
    colnames(dan) <- nvn;
    colnames(daf) <- nvf;
    dal <- cbind(daf, dan)[nv];

    dsg <- model.matrix(form, dal);	

    lss <- sapply(lsel, function(sel, dal)
              {
                  ss <- sapply(1:length(sel), function(j, sel, dal)
                           {
                               nm <- names(sel)[[j]];
                               if (is.factor(dal[, nm])) paste("'", sel[[j]], "'", sep="");
                           }, sel, dal, simplify=F);
                  names(ss) <- names(sel);
                  ss;
              }, dal, simplify=F);

    ds0 <- sapply(lss, function(ss, d)
             {
                 eval(parse(text=paste(sapply(1:length(ss), function(j, ss, d)
                        {
                            o <- ss[[j]];
                            names(o) <- names(ss)[j];
                            paste("d$", names(ss)[j], "%in%c(",
                                  paste(o, collapse=', ', sep=''), ")", sep='');
                        },  ss, dal), collapse=' & ')))
             }, d)

    if (any(apply(ds0, 1, sum) > 1)) stop("Error: Overlapping groups!");
    
    if (type=='balanced')
    {
        ds <- sapply(lss, function(ss, dal)
                 {
                     eval(parse(text=paste(sapply(1:length(ss), function(j, ss, dal)
                            {
                                o <- ss[[j]];
                                names(o) <- names(ss)[j];
                                paste("dal$", names(ss)[j], "%in%c(",
                                      paste(o, collapse=', ', sep=''), ")", sep='');
                            },  ss, dal), collapse=' & ')))
                 }, dal)

        grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
        cont <- (t(grs)%*%signs)[, 1]

    }else
    {
        dsg0 <- model.matrix(form, d);
        
        nvm <- colnames(dsg)[regexpr("\\:", colnames(dsg)) < 0];
        lnve <- sapply(ldsg, function(ds, nvm)
                   {
                       nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) == 1)];
                   }, nvm);
        lnvd <- sapply(ldsg, function(ds, nvm)
                   {
                       nvm[apply(ds[, nvm, drop=F], 2, function(o) length(unique(o)) > 1)];
                   }, nvm);
        lmed.nvd <- sapply(lnvd, function(v, dsg0) apply(dsg0[, v, drop=F], 2, mean), dsg0);

        lp <- sapply(1:length(med.nvd), function(j, nvm)
                 {
                     p <- data.frame(t(c(unique(ldsg[[j]][, lnve[[j]]]), lmed.nvd[[j]])));
                     colnames(p) <- c(lnve[[j]], lnvd[[j]]);
                     p[, nvm, drop=F];
                 }, nvm);

        grs <- t(apply(ds, 2, function(s, dsg) apply(dsg[s, , drop=F], 2, mean), dsg));
        cont <- (t(grs)%*%signs)[, 1]

    }

    if (!is.null(ngroups)) rownames(grs) <- colnames(ds0) <- ngroups;
    list(grs=grs, cont=cont, sels=ds0);

}


#################################################################################################################
#################################################################################################################
#################################################################################################################
