#########################################
## 0. Libraries and directories
#########################################
library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)
library(Rmagic)
library(Seurat)
library(parallel)
library(ggplot2)
library(ggExtra)
library(gridExtra)
library(Rmagic)
library(biomaRt)
library(RColorBrewer)
library(org.Hs.eg.db)
library(annotate)
library(dplyr)
library(grid)
library(ggridges)
library(openxlsx)
##install_version("Seurat", version = "3.2.3", repos = "http://cran.us.r-project.org")

#########################################
## 1. Import data
#########################################
mainDir <- "pathtomainDir"
load( file=paste0(mainDir, "/2.filter_mitochondrial_all/seus.RData"))

repDir <- "pathtorepDir"
tabDir <- "pathtotabDir"
figDir <- "pathtofigDir"
rouDir <- "pathtorouDir"

## In-house code - Functions 
source(paste0(rouDir, "utils.R"))
source(paste0(rouDir, "functions.r"))
library(MDA)

#########################################
## 2. Preprocess
#########################################

## get MT percentage
tit <- "allcond"
names1 <-  names(seus)
seust <- mclapply(seus[names1], function(x) {
    x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
    x
}, mc.cores=4)


## Separated by replicate
seustREP1 <- subset(seust[["all"]], features = rownames(seust[["all"]])[!grepl("^RPS|^RPL", rownames(seust[["all"]]))],
                         cells=rownames(seust[["all"]]@meta.data)[seust[["all"]]@meta.data$percent.mt < 20 &
                                                                  regexpr("rep1",rownames(seust[["all"]]@meta.data))>0])
seustREP2 <- subset(seust[["all"]], features = rownames(seust[["all"]])[!grepl("^RPS|^RPL", rownames(seust[["all"]]))],
                         cells=rownames(seust[["all"]]@meta.data)[seust[["all"]]@meta.data$percent.mt < 20 &
                                                                  regexpr("rep2",rownames(seust[["all"]]@meta.data))>0])
seustAll <- list(seust.REP1=seustREP1, seust.REP2=seustREP2)


#########################################
## 3. Normalization, dimension reduction and clustering
#########################################

## normalize by SCT
seustAll2 <- mclapply(seustAll, function(x){
  SCTransform(x, vars.to.regress = c("percent.mt"), verbose = FALSE, return.only.var.genes = FALSE)
},mc.cores=2)


## Anchors integration approach
tfeats <- SelectIntegrationFeatures(object.list = seustAll2, nfeatures = 5000)
seustAll2 <- PrepSCTIntegration(object.list = seustAll2, anchor.features = tfeats, verbose = TRUE)
anchors <- FindIntegrationAnchors(object.list = seustAll2, normalization.method = "SCT", anchor.features = tfeats, verbose = TRUE)
allfeatures <- lapply(seustAll2, row.names) %>% Reduce(intersect, .)
mi <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", verbose = TRUE, features = allfeatures)
mi <- ScaleData(object = mi, verbose = FALSE)
seust.combined <- mi

## Dimension reduction and clustering
seust.combined <- RunPCA(seust.combined, verbose = FALSE)
pdf(paste0(figDir, "elbows.pdf"))
print(ElbowPlot(seust.combined, ndims=25))
dev.off()

x <- seust.combined
x <- RunTSNE(x, dims = 1:10)
x <- RunUMAP(x, dims = 1:10)
x <- FindNeighbors(x, dims = 1:10, verbose = FALSE)
x <- FindClusters(x, verbose = FALSE, resolution=1.2)
seust.combined <-  x

## Rename samples
seust.combined$sampleid <- factor(seust.combined@meta.data$hash.ID,
                                  levels = c("parental", "MAPKi-resist","Ranolazine"))
seust.combined$sampleid_rep<- seust.combined$sampleid:as.factor(seust.combined$orig.ident)

## UMAPS representation
mycols8 <- brewer.pal(n = 10, name = "Paired")[c(9,10,3:6)]
g <- ggplotGrob(DimPlot(seust.combined, reduction="umap", group.by="sampleid_rep", dims = c(1,2),cols = mycols8) + ggtitle(tit))
ggsave(paste0(figDir, "cluster.umapall",tit,".png"),g, width=8, height=8)

g <- DimPlot(seust.combined, reduction="umap",group.by="seurat_clusters",
             cols=rainbow(length(unique(seust.combined[["seurat_clusters"]])[,1])),label = TRUE, label.size = 6)
ggsave(paste0(figDir, "cluster.umapallSHOW",tit,".png"),g, width=8, height=8)

seust.combined$readcounts <- log(apply(seust.combined$RNA[,],2,sum))

## clusters vs samples
labr <-  c("parental:rep1","parental:rep2","MAPKi-resist:rep1",
           "MAPKi-resist:rep2", "Ranolazine:rep1", "Ranolazine:rep2")

tabcomp <- function(x) sapply(labr, function(o) sum(x==o))
tabclustsall <- do.call(rbind,lapply(unique(seust.combined[["seurat_clusters"]])[,1], function(o)
    tabcomp(seust.combined[["sampleid_rep"]][,1][seust.combined[["seurat_clusters"]][,1]==o])))
rownames(tabclustsall) <- paste0("clust", unique(seust.combined[["seurat_clusters"]])[,1])
tabclustsall <- tabclustsall[paste0("clust",0:16),]

tabclustsall2 <- tabclustsall
for(i in 1:nrow(tabclustsall2))
    tabclustsall2[i,] <- tabclustsall2[i,]/sum(tabclustsall2[i,])

pdf(paste0(figDir, "cluster.perc.pdf"), width =10, height =7)
par(mar=c(5.1, 4.1, 4.1, 11.1), xpd=TRUE)
barplot(t(tabclustsall2), beside = FALSE, horiz = TRUE, las =2, col = mycols8)
legend("topright", inset=c(-0.3,0), legend=labr,
       fill=mycols8, title="Sample:rep")
dev.off()


#########################################
## 4. Impute data with magic
#########################################

seust.combined@assays$MAGIC_integrated <- magic(seust.combined, verbose=F)@assays$MAGIC_integrated

### Rambow genesets
##change path "/Volumes/biostats/acaballe/RAMBOW signatures genes.xlsx"
library("EnsDb.Hsapiens.v79")
rambowgenes <- lapply(1:6, function(k){
    aux1 <- read.xlsx("/Volumes/biostats/acaballe/RAMBOW signatures genes.xlsx",sheet=k, colNames=FALSE)
    syms <- ensembldb::select(EnsDb.Hsapiens.v79, keys= aux1[,1], keytype = "GENEID", columns = c("SYMBOL"))
    rownames(syms) <- syms[,1]
    syms[,2]
})
names(rambowgenes) <- c("NCSC","INVASION","MITFtarget", "SMC", "PIGMENTATION", "PROLIFERATIVE")


stab <- rambowgenes
allgenes <- unique(unlist(stab))

## MAGIC INDIVIDUAL GENS
p <- seust.combined
for(i in allgenes){
        if(i %in% rownames(p@assays$MAGIC_integrated)) {
            z <- as.numeric(p@assays$MAGIC_integrated[i,])
            p[[paste0(i, "_mag")]] <- z
            p[[paste0(i, "_mag_scale")]] <- scale(z)
        } else {
            p[[paste0(i, "_mag")]] <- rep(0, ncol(p))
        }
    }
seust.combined <-   p

## MAGIC SIGNATURES
cmpop <- lapply(stab, function(x) colMeans(as.matrix(seust.combined@assays$integrated[rownames(seust.combined@assays$integrated) %in% x,])))
for(j in names(cmpop)) seust.combined[[j]] <- cmpop[[j]]
cmpop.mag <- lapply(stab, function(x) colMeans(as.matrix(seust.combined@assays$MAGIC_integrated[rownames(seust.combined@assays$MAGIC_integrated) %in% x,])))

for(j in names(cmpop.mag)) seust.combined[[paste0(j, "_mag")]] <- cmpop.mag[[j]]
for(j in names(cmpop.mag)) seust.combined[[paste0(j, "_mag_scale")]] <- scale(cmpop.mag[[j]])
sel <- names(stab)
sel.mag <- paste0(names(stab), "_mag")
sel.mag.scale <- paste0(names(stab), "_mag_scale")

## Plots collapsed (umaps)
K <- 1
cols <- c(colorRamps::matlab.like2(20)[1:13], "deeppink2", "deeppink3", "deeppink4")
g <- lapply(sel.mag[1:4], function(i) FeaturePlot(seust.combined, combine=TRUE, reduction="umap", features=i,
                                                        pt.size=0.7, ncol=4, cols=cols, label = TRUE)+ ggtitle(paste0(i)))
g[[5]] <- DimPlot(seust.combined, reduction="umap", group.by="seurat_clusters", cols=rainbow(length(unique(seust.combined[["seurat_clusters"]])[,1])), label = TRUE)
g[[6]] <- DimPlot(seust.combined, reduction="umap", group.by="sampleid_rep", cols=mycols8, label = FALSE)
ggsave(paste0(figDir, "signatures.umap_all1",tit,".png"), marrangeGrob(grobs = g, layout_matrix = t(rbind(c(1,2),c(3,4),c(5,6)))), width=25, height=12, dpi = 120)


cols <- c(colorRamps::matlab.like2(20)[1:13], "deeppink2", "deeppink3", "deeppink4")
g <- lapply(sel.mag[2:6], function(i) FeaturePlot(seust.combined, combine=TRUE, reduction="umap", features=i,
                                                        pt.size=0.7, ncol=4, cols=cols, label = TRUE)+ ggtitle(paste0(i)))
g[[5]] <- DimPlot(seust.combined, reduction="umap", group.by="seurat_clusters", cols=rainbow(length(unique(seust.combined[["seurat_clusters"]])[,1])), label = TRUE)
g[[6]] <- DimPlot(seust.combined, reduction="umap", group.by="sampleid_rep", cols=mycols8, label = FALSE)
ggsave(paste0(figDir, "signatures.umap_all2",tit,".png"), marrangeGrob(grobs = g, layout_matrix = t(rbind(c(1,2),c(3,4),c(5,6)))), width=25, height=12, dpi = 120)


## Plots collapsed (violins)
g <- lapply(sel.mag, function(k)  ggplotGrob(VlnPlot(seust.combined, features =k, pt.size = 0.0001, group.by = "sampleid_rep", cols = mycols8)
                + ggtitle(paste0("Violin: expression ","_",k))))
ggsave(paste0(figDir, "sample_rep_signatures.png"), marrangeGrob(grobs = g, layout_matrix =  rbind(c(1,2,3,4),c(5,6,7,8))), width=24, height=8)


for (o in c(names(stab))){
    dir.create(resp <- paste0(figDir, o,"plots_signatures/"))
    for(i in stab[[o]]){
        i <- sub("-",".",i, fixed = TRUE)
        data1 <-  data.frame(seust.combined[[paste0(i, "_mag")]],seust.combined$sampleid_rep,seust.combined$seurat_clusters )
        colnames(data1) <- c("y","sample", "clusters")

        g<-list()
        g[[1]] <- ggplotGrob(FeaturePlot(seust.combined, combine=TRUE, reduction="umap", features=paste0(i,"_mag"),
                                         pt.size=0.3, ncol=4, cols=cols, label = TRUE)+ ggtitle(paste0(i)))
        g[[2]] <- DimPlot(seust.combined, reduction="umap", group.by="seurat_clusters",
                          cols=rainbow(length(unique(seust.combined[["seurat_clusters"]])[,1])), label = TRUE)
        g[[3]] <- DimPlot(seust.combined, reduction="umap", group.by="sampleid_rep",
                          cols=mycols8, label = FALSE)

        g[[4]]  <- ggplot(
            data1,
            aes(x =y, y =sample, fill = sample))+
                geom_density_ridges(
                    aes(point_color = sample, point_fill = sample, point_shape = sample),
                    alpha = .2, point_shape = "|", point_size = 1,  point_alpha = 1, jittered_points = TRUE) +
                        scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23,24)) +
                            labs(title = paste0(""))+ xlab(i) + theme(legend.position = "none") + scale_fill_manual(values = mycols8) +
                                scale_color_manual(values = mycols8) +   scale_discrete_manual("point_color", values = mycols8, guide = "none")
        ggsave(paste0(resp, "clusters_signatures_and_violin_",i,tit,"all.png"), marrangeGrob(grobs = g, layout_matrix = rbind(c(1,2),c(3,4))),
               width=12, height=10, dpi=100)

    }
}

##  Plots individual other interesting genes
sel.mag <- c("CD36", "DUSP6","MTAP","SLC38A2","MAT2A","CD274","MITF","AXL","EDNRB", "DCT", "RAB27A", "PAX3", "SOX10", "TYR", "TRMP1","TYRP1")
sel.mag <- sel.mag[sel.mag%in%rownames(seust.combined)]

for(k in sel.mag){
    seust.combined[[paste0(k, "_mag")]] <- as.numeric(seust.combined@assays$MAGIC_integrated[k,])
}

K <- 1
cols <- c(colorRamps::matlab.like2(20)[1:18], "deeppink2", "deeppink3", "deeppink4")
dir.create(resp <- paste0(figDir, "genes_ind","plots/"))

for(i in sel.mag){
    data1 <-  data.frame(seust.combined[[paste0(i, "_mag")]],seust.combined$sampleid_rep,seust.combined$seurat_clusters )
    colnames(data1) <- c("y","sample", "clusters")

    g<-list()
    g[[1]] <- ggplotGrob(FeaturePlot(seust.combined, combine=TRUE, reduction="umap", features=paste0(i,"_mag"),
                                     pt.size=0.7, ncol=4, cols=cols, label = TRUE)+ ggtitle(paste0(i)))
    g[[2]] <- DimPlot(seust.combined, reduction="umap", group.by="seurat_clusters",
                      cols=rainbow(length(unique(seust.combined[["seurat_clusters"]])[,1])), label = TRUE)
    g[[3]] <- DimPlot(seust.combined, reduction="umap", group.by="sampleid_rep",
                  cols=mycols8, label = FALSE)

    g[[4]]  <- ggplot(
  data1,
  aes(x =y, y =sample, fill = sample))+
      geom_density_ridges(
    aes(point_color = sample, point_fill = sample, point_shape = sample),
    alpha = .2, point_shape = "|", point_size = 1,  point_alpha = 1, jittered_points = TRUE) +
 scale_discrete_manual(aesthetics = "point_shape", values = c(21, 22, 23,24)) +
     labs(title = paste0(""))+ xlab(i) + theme(legend.position = "none") + scale_fill_manual(values = mycols8) +
         scale_color_manual(values = mycols8) +   scale_discrete_manual("point_color", values = mycols8, guide = "none")
    ggsave(paste0(resp, "signatures.umap_all",i,tit,".png"), marrangeGrob(grobs = g, layout_matrix = rbind(c(1,2),c(3,4))), width=15, height=12, dpi = 60)
}


###################################################
########### 5. MARKERS BY CLUSTERS ################
###################################################

## Markers by clusters
pbmc.markers <-  FindAllMarkers(seust.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers, sep ="\t", quote = FALSE, row.names = FALSE, file = paste0(tabDir, "clust.markers",tit,".xls"))

sel <- 	top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
x <- as.data.frame(top10)
x <- cbind(gene = x[,7], x[,2:4], format(x[,c(1,5)],digits = 3), cluster = x[,"cluster"])


## Biological enrichment (rambow genesets)
gsets <- rambowgenes
GENES <- rownames(seust.combined@assays$integrated)
forD <- lapply(levels(pbmc.markers$cluster), function(s) pbmc.markers[pbmc.markers$cluster == s,"gene"])
names(forD) <- paste0("cluster ",levels(pbmc.markers$cluster))
GENESETS2 <- list(gsets)
gsets <- "signatures_interest"
names(GENESETS2) <- gsets

dir.create(res <- paste0(tabDir, "GOEnrich_rambow/"))
overlaps <- runHyper(GENESETS2, gsets, forD, res)
overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)

## Biological enrichment (GO/KEGG terms)

### change path "/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/" to v3 hallmarks, GO terms! 
gsets <- paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/",
               list.files(paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/", sep='')), sep='');
gsets <- gsets[regexpr("Clus.", gsets) < 0];
gsets <- gsets[regexpr("expevidence.", gsets) < 0];
gsets <- gsets[regexpr("differen.", gsets) < 0];
gsets <- gsets[regexpr("v3.", gsets) > 0];

GENES <- rownames(seust.combined@assays$integrated)
forD <- lapply(levels(pbmc.markers$cluster), function(s) pbmc.markers[pbmc.markers$cluster == s,"gene"])
names(forD) <- paste0("cluster ",levels(pbmc.markers$cluster))
GENESETS2 <- prepareGS(gsets, GENES)

dir.create(res <- paste0(tabDir, "GOEnrich_GO/"))
overlaps <- runHyper(GENESETS2, gsets, forD, res)
overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)



###################################################
########### 6. MARKERS BY SAMPLE   ################
###################################################
seust.combined2 <- seust.combined
seust.combined2[["seurat_clusters"]] <- seust.combined2$sampleid
Idents(object = seust.combined2) <- seust.combined2$sampleid

pbmc.markers2 <-  FindAllMarkers(seust.combined2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(pbmc.markers2, sep ="\t", quote = FALSE, row.names = FALSE, file = paste0(tabDir, "clust.markers",tit,"_sampleid.xls"))

sel <- 	top10 <- pbmc.markers2 %>% group_by(cluster) %>% top_n(n = 20, wt = avg_logFC)
x <- as.data.frame(top10)
x <- cbind(gene = x[,7], x[,2:4], format(x[,c(1,5)],digits = 3), cluster = x[,"cluster"])


## Biological enrichment (rambow genesets)
gsets <- rambowgenes
GENES <- rownames(seust.combined@assays$integrated)
forD <- lapply(levels(pbmc.markers2$cluster), function(s) pbmc.markers2[pbmc.markers2$cluster == s,"gene"])
names(forD) <- paste0("cluster ",levels(pbmc.markers2$cluster))
GENESETS2 <- list(gsets)
gsets <- "signatures_interest"
names(GENESETS2) <- gsets

dir.create(res <- paste0(tabDir, "GOEnrich_samples_rambow/"))
overlaps <- runHyper(GENESETS2, gsets, forD, res)
overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)

## Biological enrichment (GO/KEGG terms)

### change path "/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/" to v3 hallmarks, GO terms! 
gsets <- paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/",
               list.files(paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/", sep='')), sep='');
gsets <- gsets[regexpr("Clus.", gsets) < 0];
gsets <- gsets[regexpr("expevidence.", gsets) < 0];
gsets <- gsets[regexpr("differen.", gsets) < 0];
gsets <- gsets[regexpr("v3.", gsets) > 0];

GENES <- rownames(seust.combined@assays$integrated)
forD <- lapply(levels(pbmc.markers2$cluster), function(s) pbmc.markers2[pbmc.markers2$cluster == s,"gene"])
names(forD) <- paste0("cluster ",levels(pbmc.markers2$cluster))
GENESETS2 <- prepareGS(gsets, GENES)

dir.create(res <- paste0(tabDir, "GOEnrich_samples_GO/"))
overlaps <- runHyper(GENESETS2, gsets, forD, res)
overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)



################################################################
######## 7. Rambow top cells vs other pathways #################
################################################################
library(reshape2)

## discretization of rambow expression (Selection of top 10% most expressed cells)
stab <- rambowgenes
y <- seust.combined
cmpop <-  lapply(stab, function(x) colMeans(as.matrix(y@assays$integrated[rownames(y@assays$integrated) %in% x,])))
for(j in names(cmpop)) y[[j]] <- cmpop[[j]]
cmpop.mag <-  lapply(stab, function(x) colMeans(as.matrix(y@assays$MAGIC_integrated[rownames(y@assays$MAGIC_integrated) %in% x,])))
for(j in names(cmpop.mag)) y[[paste0(j, "_mag")]] <- cmpop.mag[[j]]

for(j in names(cmpop.mag)){
    y[[paste0(j, "_mag_discrete")]] <- ifelse(cmpop.mag[[j]]> quantile(cmpop.mag[[j]],0.9),"high","low")
}

## Visualization of discretized signature expression
K <- 1
cols <- c(colorRamps::matlab.like2(20)[1:13], "deeppink2", "deeppink3", "deeppink4")
for(i in names(stab)){
    g <- list()
    g[["Magic"]] <-  FeaturePlot(y, combine=TRUE, reduction="umap", features=paste0(i,"_mag"),
                      pt.size=0.7, ncol=4, cols=cols, label = FALSE )+ ggtitle(paste0(i))
    g[["Discrete"]] <-  DimPlot(y, combine=TRUE, reduction="umap", group.by=paste0(i, "_mag_discrete"),
                                pt.size=0.7, ncol=4, cols=c("gray70","blue"), label = FALSE,
                                order=c("high","low"))+ ggtitle("")
    g[["violin"]] <- VlnPlot(y, features =paste0(i,"_mag"), pt.size = 0, group.by = "sampleid_rep", cols = mycols8) + ggtitle("")

    df <- table(data.frame(y[["sampleid_rep"]], y[[paste0(i, "_mag_discrete")]][,1]))
    df <- data.frame(ncells = c(df[,1],df[,2]) /c(rep(df[,1]+df[,2],2)), sample  = factor(rep(rownames(df),2),levels= rownames(df)),
                     state = rep(colnames(df),each =6))
    g[["Samples"]] <- ggplot(df, aes(x= sample, y = ncells, fill = state,label = scales::percent(ncells))) + geom_bar(stat="identity") +theme_bw() +
        coord_flip() +xlab("percent cells") + scale_fill_manual(values=c("blue","gray70"))


    ggsave(paste0(figDir, "signaturesrambow.umap_",i,tit,"discrete.pdf"),
           marrangeGrob(grobs = g, layout_matrix = t(rbind(c(1,2),c(3,4)))), width=12, height=10, dpi = 120)
}


### Markers and enrichment

### change path "/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/" to v3 hallmarks, GO terms! 
gsets <- paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/",
               list.files(paste("/Volumes/biostats/databases/BroadGSEA/genesets/hsapiens/", sep='')), sep='');
gsets <- gsets[regexpr("Clus.", gsets) < 0];
gsets <- gsets[regexpr("expevidence.", gsets) < 0];
gsets <- gsets[regexpr("differen.", gsets) < 0];
gsets <- gsets[regexpr("v3.", gsets) > 0];

res <- paste0(tabDir, "GOEnrich_corr2rambow/")
system(paste0("mkdir ", res))

###### pigmentation
i <- names(stab)[2]
i2 <- "pigmentation"

for(S in c("MAPKi-resist","Ranolazine")){
    yp <- subset(y, features = rownames(y),
                    cells=rownames(y@meta.data)[y$sampleid==S])
    markers <- FindMarkers(yp, ident.1 = "high", group.by = paste0(i, "_mag_discrete"), logfc.threshold=0.2, min.pct=0.2, only.pos=TRUE)
    markersG <- rownames(markers)[1:min(200,nrow(markers))]
    c1 <- cor(t(yp@assays$MAGIC_integrated@data), yp[[paste0(i, "_mag")]])

    ## enrichment
    GENES <- rownames(yp@assays$integrated)
    forD <- list(markers_discrete = markersG, markers_continuous = names(head(c1[order(-c1),],200)))
    names(forD) <- paste0(names(forD),"_",i2)

    GENESETS2 <- prepareGS(gsets, GENES)
    dir.create(res <- paste0(figDir,i2, sub("-","_",S),"/"))
    overlaps <- runHyper(GENESETS2, gsets, forD, res)
    overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)
    htmlobj <- buildResHTML(res, forD)
    colnames(htmlobj$xout) <- names(forD)
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,".html"),
                   links=htmlobj$links, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
    links2 <- htmlobj$links
    links2[[1]] <- sub("xls","html",links2[[1]])
    links2[[2]] <- sub("xls","html",links2[[2]])
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,"_HTML.html"),
                   links=links2, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)

}

###### mitosis
i <- names(stab)[1]
i2 <- "mitosis"

for(S in c("MAPKi-resist","Ranolazine","parental")){
    yp <- subset(y, features = rownames(y),
                    cells=rownames(y@meta.data)[y$sampleid==S])
    markers <- FindMarkers(yp, ident.1 = "high", group.by = paste0(i, "_mag_discrete"), logfc.threshold=0.2, min.pct=0.2, only.pos=TRUE)
    markersG <- rownames(markers)[1:min(200,nrow(markers))]
    c1 <- cor(t(yp@assays$MAGIC_integrated@data), yp[[paste0(i, "_mag")]])

    ## enrichment
    GENES <- rownames(yp@assays$integrated)
    forD <- list(markers_discrete = markersG, markers_continuous = names(head(c1[order(-c1),],200)))
    names(forD) <- paste0(names(forD),"_",i2)

    GENESETS2 <- prepareGS(gsets, GENES)
    dir.create(res <- paste0(figDir,i2, sub("-","_",S),"/"))
    overlaps <- runHyper(GENESETS2, gsets, forD, res)
    overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)
    htmlobj <- buildResHTML(res, forD)
    colnames(htmlobj$xout) <- names(forD)
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,".html"),
                   links=htmlobj$links, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
    links2 <- htmlobj$links
    links2[[1]] <- sub("xls","html",links2[[1]])
    links2[[2]] <- sub("xls","html",links2[[2]])
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,"_HTML.html"),
                   links=links2, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
}

###### invasive
i <- names(stab)[3]
i2 <- "invasive"

for(S in c("MAPKi-resist","Ranolazine")){
    yp <- subset(y, features = rownames(y),
                    cells=rownames(y@meta.data)[y$sampleid==S])
    markers <- FindMarkers(yp, ident.1 = "high", group.by = paste0(i, "_mag_discrete"), logfc.threshold=0.2, min.pct=0.2, only.pos=TRUE)
    markersG <- rownames(markers)[1:min(200,nrow(markers))]
    c1 <- cor(t(yp@assays$MAGIC_integrated@data), yp[[paste0(i, "_mag")]])

    ## enrichment
    GENES <- rownames(yp@assays$integrated)
    forD <- list(markers_discrete = markersG, markers_continuous = names(head(c1[order(-c1),],200)))
    names(forD) <- paste0(names(forD),"_",i2)

    GENESETS2 <- prepareGS(gsets, GENES)
    dir.create(res <- paste0(figDir,i2, sub("-","_",S),"/"))
    overlaps <- runHyper(GENESETS2, gsets, forD, res)
    overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)
    htmlobj <- buildResHTML(res, forD)
    colnames(htmlobj$xout) <- names(forD)
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,".html"),
                   links=htmlobj$links, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
    links2 <- htmlobj$links
    links2[[1]] <- sub("xls","html",links2[[1]])
    links2[[2]] <- sub("xls","html",links2[[2]])
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,"_HTML.html"),
                   links=links2, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)

}

###### NSCS
i <- names(stab)[4]
i2 <- "NSCS"

for(S in c("MAPKi-resist")){
    yp <- subset(y, features = rownames(y),
                    cells=rownames(y@meta.data)[y$sampleid==S])
    markers <- FindMarkers(yp, ident.1 = "high", group.by = paste0(i, "_mag_discrete"), logfc.threshold=0.2, min.pct=0.2, only.pos=TRUE)
    markersG <- rownames(markers)[1:min(200,nrow(markers))]
    c1 <- cor(t(yp@assays$MAGIC_integrated@data), yp[[paste0(i, "_mag")]])

    ## enrichment
    GENES <- rownames(yp@assays$integrated)
    forD <- list(markers_discrete = markersG, markers_continuous = names(head(c1[order(-c1),],200)))
    names(forD) <- paste0(names(forD),"_",i2)

    GENESETS2 <- prepareGS(gsets, GENES)
    dir.create(res <- paste0(figDir,i2, sub("-","_",S),"/"))
    overlaps <- runHyper(GENESETS2, gsets, forD, res)
    overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)
    htmlobj <- buildResHTML(res, forD)
    colnames(htmlobj$xout) <- names(forD)
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,".html"),
                   links=htmlobj$links, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
    links2 <- htmlobj$links
    links2[[1]] <- sub("xls","html",links2[[1]])
    links2[[2]] <- sub("xls","html",links2[[2]])
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,"_HTML.html"),
                   links=links2, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)

}


###### MITF.targets.signature
i <- names(stab)[5]
i2 <- "MITF.targets.signature"

for(S in c("MAPKi-resist","Ranolazine")){
    yp <- subset(y, features = rownames(y),
                    cells=rownames(y@meta.data)[y$sampleid==S])
    markers <- FindMarkers(yp, ident.1 = "high", group.by = paste0(i, "_mag_discrete"), logfc.threshold=0.2, min.pct=0.2, only.pos=TRUE)
    markersG <- rownames(markers)[1:min(200,nrow(markers))]
    c1 <- cor(t(yp@assays$MAGIC_integrated@data), yp[[paste0(i, "_mag")]])

    ## enrichment
    GENES <- rownames(yp@assays$integrated)
    forD <- list(markers_discrete = markersG, markers_continuous = names(head(c1[order(-c1),],200)))
    names(forD) <- paste0(names(forD),"_",i2)

    GENESETS2 <- prepareGS(gsets, GENES)
    dir.create(res <- paste0(figDir,i2, sub("-","_",S),"/"))
    overlaps <- runHyper(GENESETS2, gsets, forD, res)
    overlaps2 <- runHyperhtml(GENESETS2, gsets, forD, res)
    htmlobj <- buildResHTML(res, forD)
    colnames(htmlobj$xout) <- names(forD)
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,".html"),
                   links=htmlobj$links, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
    links2 <- htmlobj$links
    links2[[1]] <- sub("xls","html",links2[[1]])
    links2[[2]] <- sub("xls","html",links2[[2]])
    write.html.mod(htmlobj$xout, file=paste0(res, "summary",i2,S,"_HTML.html"),
                   links=links2, tiny.pic=vector(mode="list", length=ncol(htmlobj$xout)),title=tit)
}


