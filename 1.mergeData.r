#########################################
## 0. Libraries and directories
#########################################
library(Matrix)
library(parallel)
library(hwriter)
library(Seurat)
library(parallel)

mainDir <- "path2maindir"
dataDir <- paste0(mainDir, "data/")
repDir <- paste0(mainDir, "reports/1.mergeData/")
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

load(paste0(mainDir, "reports/1.readData_rep2/seus_r2.RData"))
seus2 <- seus
load(paste0(mainDir, "reports/1.readData_rep1/seus_r1.RData"))

for(i in names(seus))
    seus[[i]] <- merge(seus[[i]], y = seus2[[i]], add.cell.ids = c("rep1", "rep2"), project = "Melanoma")

head(colnames(seus[[i]]))
seus[[i]]$orig.ident <- as.factor(substr(colnames(seus[[i]]),1,4))
save(seus, file = paste0(repDir, "seus.RData"))
