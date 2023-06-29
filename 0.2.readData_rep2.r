#########################################
## 0. Libraries and directories
#########################################
library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)
library(parallel)

mainDir <- "path2maindir/"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/1.readData_rep2/")
figDir <- paste0(repDir, "figs/")
tabDir <- paste0(repDir, "tables/")

dir.create(repDir)
dir.create(figDir)
dir.create(tabDir)

mycols <- c(1, '#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506')

corner <- function(x) x[1:5,1:5]


#########################################
## 1. Import data
#########################################

#### Read matrix format

barcode.path <- paste0(dataDir, "seq_rep2/barcodes.tsv")
features.path <- paste0(dataDir, "seq_rep2/features.tsv")
matrix.path <- paste0(dataDir, "seq_rep2/matrix.mtx")

mats <- mclapply(1:length(barcode.path), function(i){
    mat <- readMM(file = matrix.path[i])
    gene.names = read.delim(features.path[i],
        header = FALSE,
        stringsAsFactors = FALSE)
    barcode.names = read.delim(barcode.path[i],
        header = FALSE,
        stringsAsFactors = FALSE)
    colnames(mat) = barcode.names$V1
    rownames(mat) = gene.names$V1
    mat
}, mc.cores=1)


### Check UMI filters
sel <- 36602:36605
sampleiddata <-  (mats[[1]][sel,])
rownames(sampleiddata) = c("parental","MAPKi_resist","Ranolazine","Etomoxir")
table(apply(sampleiddata,1,which.max))

filters <- list(1000, 20000, 100)
names(filters) <- c("minUMI",  "maxGene", "minGene")

ncmus <- mclapply(mats, function(mat){
    cs <- Matrix::colSums(mat[-sel,])
    cg <- Matrix::colSums(mat[-sel,] > 0)
    sapply(seq(1000, to=3000, by=100), function(mu) sum(cg < filters$maxGene & cg > filters$minGene & cs > mu))
}, mc.cores=1)
ncm <- do.call(rbind, ncmus)
colnames(ncm) <- seq(1000, to=3000, by=100)

nmax <- max(unlist(ncmus))
nmin <- min(unlist(ncmus))

library(Seurat)

### Set umi filter at 1000
sel2 <- 36602:36605

mat <- mats[[1]]
cs <- Matrix::colSums(mat[-sel2,])
cg <- Matrix::colSums(mat[-sel2,] > 0)
sel <- cs > filters$minUMI & cg < filters$maxGene & cg > filters$minGene
rs <- Matrix::rowSums(mat[-sel2,])
matf <- mat[which(rs > 0), sel]
sampleiddata2 <- sampleiddata[,sel]
matf <- list(matf=matf, stats=list(cs=cs, cg=cg), sampleiddata = sampleiddata2)


mats <- list(matf[[1]])
sampleiddata <- matf$sampleiddata

### Name genes by symbol

library("EnsDb.Hsapiens.v79")

syms <- ensembldb::select(EnsDb.Hsapiens.v79, keys= rownames(mats[[1]]), keytype = "GENEID", columns = c("SYMBOL"))
syms <- syms[match(unique(syms[,2]), syms[,2]),]
mats[[1]] <- mats[[1]][match(syms[,1], rownames(mats[[1]])),]
rownames(mats[[1]]) <- syms[,2]


### Create Seurat object and add hash info

seus <- mclapply(mats, function(x){
    seu <- CreateSeuratObject(x)
    seu[["hash"]] <-  CreateAssayObject(counts = sampleiddata)
    seu <- NormalizeData(seu, assay = "hash", normalization.method = "CLR")
    seu <- HTODemux(seu, assay = "hash", positive.quantile = 0.99)
    seu
}, mc.cores=1)
names(seus) <- "all"

pdf(paste0(figDir, "maxID.ridge_r2.pdf"), width=25, height=5)
for (i in names(seus)) {
    seu <- seus[[i]]
    Idents(seu) <- "hash_maxID"
    print(RidgePlot(seu, assay = "hash", features = rownames(seu[["hash"]])[1:4], ncol = 4))
}
dev.off()

pdf(paste0(figDir, "hash.ID.ridge_r2.pdf"), width=25, height=8)
for (i in names(seus)) {
    seu <- seus[[i]]
    Idents(seu) <- "hash.ID"
    print(RidgePlot(seu, assay = "hash", features = rownames(seu[["hash"]])[1:4], ncol = 4))
}
dev.off()

pdf(paste0(figDir, "hash_classification.ridge_r2.pdf"), width=55, height=8)
for (i in names(seus)) {
    seu <- seus[[i]]
    Idents(seu) <- "hash_classification"
    print(RidgePlot(seu, assay = "hash", features = rownames(seu[["hash"]])[1:4], ncol = 4))
}
dev.off()

seu <- seus[[1]]
Idents(seu) <- "hash_classification.global"
seu <- subset(seu, idents = "Negative", invert = TRUE)
seu <- subset(seu, idents = "Doublet", invert = TRUE)

seus <- mclapply(seus, function(seu) {
    Idents(seu) <- "hash_classification.global"
    seu <- subset(seu, idents = "Negative", invert = TRUE)
    seu <- subset(seu, idents = "Doublet", invert = TRUE)
    SplitObject(seu, split.by = "hash_classification")
}, mc.cores=1)

seus <- unlist(seus, recursive=FALSE)
names(seus) <- sub("all.","",names(seus))
seus[["all"]] <- seu

save(seus, file=paste0(repDir, "seus_r2.RData"))
