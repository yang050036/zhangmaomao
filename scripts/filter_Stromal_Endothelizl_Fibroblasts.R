library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)
library(metap)
library(SingleR)
library(egg)
library(reshape2)
library(gtools)
mus <- readRDS("ImmGenData.rds")
filters <- c('Stromal cells','Endothelial cells','Fibroblasts')
immune.combine$filter <- unlist(lapply(immune.combine$SingleR.cluster.lables,function(x){x %in% filters}))
split <- SplitObject(immune.combine,split.by = 'filter')
filter <- split$'FALSE'

DefaultAssay(filter) <- "RNA"
all.genes <- row.names(filter)
DefaultAssay(filter) <- "integrated"
filter <- ScaleData(filter,features = all.genes,verbose = FALSE)
filter <- RunPCA(filter,npcs = 30,verbose = FALSE)
filter <- FindNeighbors(filter,reduction = 'pca',dims = 1:20)
filter <- FindClusters(filter,resolution = 0.5)
filter <- RunTSNE(filter,reduction = "pca",dims = 1:20)
filter <- RunUMAP(filter,reduction = "pca",dims = 1:20)
rds.ex <- as.SingleCellExperiment(rds)
singlr = SingleR(test=rds.ex,assay.type.test=1,ref=mus,labels=mus$label.main,method = "cluster",clusters = rds.ex$seurat_clusters)
rds[["celltype"]] <- 
  singlr$labels[match(rds.ex[["seurat_clusters"]],row.names(singlr))]

filter.rename <- RenameIdents(filter,'0'='M_1','1'='M_2','2'='Monocytes','3'='B',
                              '4'='N_1','5'='M_3','6'='N_2','7'='M_4','8'='M_5',
                              '9'='M_6','10'='T','11'='M_7','12'='M_8','13'='M_9','14'='NK',
                              '15'='DC_1','16'='DC_2')