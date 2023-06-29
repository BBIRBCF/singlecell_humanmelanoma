#########################################
## 0. Libraries and directories
#########################################

library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)
library(numbers)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/2.filter_mitochondrial_all/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')
mycols19 <- c('#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', rainbow(10))
corner <- function(x) x[1:5,1:5]

#########################################
## 1. Import data
#########################################

## Filtered matrices

load(paste0(mainDir, "reports/1.mergeData/seus.RData"))


#########################################
## 3. Normalize samples
#########################################

seus <- lapply(seus, function(sot) {
    sot[["percent.mt"]] <- PercentageFeatureSet(sot, pattern = "^MT-")
    sot[["RPL"]] <- PercentageFeatureSet(sot, pattern = "^RPL")
    sot[["RPS"]] <- PercentageFeatureSet(sot, pattern = "^RPS")
    sot[["RPA"]] <- PercentageFeatureSet(sot, pattern = "^RPS|^RPL")
    sot
})
qcs <- lapply(c("percent.mt", "RPS", "RPL", "RPA"), function(x) {
    y <- do.call(rbind,lapply(names(seus), function(i) data.frame(seus[[i]][[x]][,1], i,seus[[i]][["orig.ident"]][,1])))
    colnames(y) <- c("value", "sample", "rep")
    y
})
names(qcs) <- c("percent.mt", "RPS", "RPL", "RPA")


mitopc <- sapply(seq(2,70, by=2), function(r){
    tt <- do.call(rbind, tapply(qcs[[1]][,1], qcs[[1]][,2]:qcs[[1]][,3], function(x) table(x < r)))
    (tt/rowSums(tt))[,2]
})
colnames(mitopc) <- seq(2,70, by=2)
mitopc[mitopc==0.5] <- NA

tt <- do.call(rbind, tapply(qcs[[1]][,1], qcs[[1]][,2]:qcs[[1]][,3], function(x) table(x < 20)))
sink(file=paste0(tabDir, "mito_perc_20pc.txt"))
tt
tt/rowSums(tt)
sink()

library(ggplot2)
pdf(paste0(figDir, "boxplot.mt.rpa.pdf"), width=16, height=12)
par(mfrow=c(2,2))
for(i in c("percent.mt", "RPS", "RPL", "RPA")) {
    qcs[[i]]$sample_rep <- qcs[[i]]$sample:qcs[[i]]$rep
    g <- ggplot(aes(x=sample_rep, y=value, colour=sample_rep), data=qcs[[i]])
    print(g + geom_boxplot(outlier.shape = NA) + ggtitle(i))
}
dev.off()

ncolini <- sapply(seus, ncol)
save(ncolini, file=paste0(repDir, "ncolini.RData"))

### Remove RPA and filter by mito perc

seusf <- mclapply(seus, function(s) subset(s, features = rownames(s)[!grepl("^RPS|^RPL", rownames(s))],
                  cells=rownames(s@meta.data)[s@meta.data$percent.mt < 20]), mc.cores=6)

nco <- sapply(seus, function(x) table(x[["orig.ident"]][,1]))
ncf <- sapply(seusf, function(x) table(x[["orig.ident"]][,1]))

pdf(paste0(figDir, "mito_perc_NC.pdf"))
plot(nco, ncf, pch=rep(c(21,22),5), bg=mycols6)
abline(0,1,lty=2)
legend("topleft", legend=rownames(mitopc), pt.bg=mycols6, pch = rep(c(21,22),5), cex=0.8)
dev.off()

ncolitof <- sapply(seusf, ncol)
save(ncolitof, file=paste0(repDir, "ncolitof.RData"))

df <- data.frame(namesID = rep(colnames(ncf), each = 2), rep = rep(rownames(ncf),5),inicells = as.numeric(nco), mitofilt = as.numeric(ncf))
write.table(df, file = paste0(repDir, "tables/mitofilt.xls"))

seus <- seusf
save(seus, file=paste0(repDir, "seus.RData"))


