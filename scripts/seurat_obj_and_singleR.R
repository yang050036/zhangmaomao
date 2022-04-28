#!/usr/bin/Rscript

library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)
library(metap)
library(SingleR)
library(egg)
library(reshape2)
library(RColorBrewer)
library(SingleCellExperiment)
library(svglite)
library(pheatmap)
library(monocle3)
library(scater)
#合并两个样本的数据
#人类免疫细胞图谱：DatabaseImmuneCellExpressionData
#人类原初细胞图谱：HumanPrimaryCellAtlasData
#小鼠批量表达数据：ImmGenData
#人类免疫细胞的RNA-seq数据：MonacoImmuneData
#已分类细胞群体的小鼠批量表达数据：MouseRNAseqData
#人类血细胞微列阵表达数据：NovershternHematopoieticData

ggplot2_color <- c()
mus <- readRDS("/mnt/beegfs/Research/Database/Seurat/ImmGenData.rds")
setwd('D:\\project\\scRNA\\zhangmaomao')
files <- list.files('./new/',pattern = 'd')
mixedsort(files)
filelist <- lapply(files,function(file){
  CreateSeuratObject(counts = Read10X(file.path("new/", file)),
                     project = file,min.cells = 3,min.features = 200)
})

day <- merge(filelist[[1]],y = c(filelist[[2]],filelist[[3]],filelist[[4]],filelist[[5]],filelist[[6]]),add.cell.ids = files,project = 'day')
#day <- merge(filelist[[1]],y = c(filelist[[2]]),add.cell.ids = c(files[1],files[2]),project = 'day')

day[['percent.mt']] <- PercentageFeatureSet(day,pattern = "^mt-")
Vp <- VlnPlot(day,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
plot1 <- FeatureScatter(day,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(day,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
point <- plot1 + plot2
day <- subset(day,subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt<50)
day.list <- SplitObject(day,split.by = 'orig.ident')
day.list <- lapply(X = day.list,FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method = 'vst',nfeatures = 2000)
})

#寻找两个样本的锚点
immune.anchors <- FindIntegrationAnchors(object.list = day.list,dims = 1:20)
#整合两个数据
immune.combine <- IntegrateData(anchorset = immune.anchors,dims = 1:20)



#immune.combine=SingleR(test=immune.combine,assay.type.test=1,ref=mus,labels=mus$label.main)

#和单样本一样的常规分析
all.genes <- row.names(immune.combine)
DefaultAssay(immune.combine) <- "integrated"
immune.combine <- ScaleData(immune.combine,features = all.genes,verbose = FALSE)
immune.combine <- RunPCA(immune.combine,npcs = 30,verbose = FALSE)
immune.combine <- RunTSNE(immune.combine,reduction = "pca",dims = 1:20)
immune.combine <- RunUMAP(immune.combine,reduction = "pca",dims = 1:20)
immune.combine <- FindNeighbors(immune.combine,reduction = 'pca',dims = 1:20)
immune.combine <- FindClusters(immune.combine,resolution = 0.5)
ep <- as.SingleCellExperiment(immune.combine)
singlr = SingleR(test=ep,assay.type.test=1,ref=mus,labels=mus$label.main,method = "cluster",clusters = ep$seurat_clusters)
immune.combine[["SingleR.cluster.lables"]] <- 
  singlr$labels[match(immune.combine[[]][["seurat_clusters"]],row.names(singlr))]

#for(i in 1:length(levels(immune.combine$seurat_clusters))){
immune.rename1 <- RenameIdents(immune.combine,'0'='Macrophages','1'='Fibroblasts','2'='Fibroblasts','3'='Fibroblasts',
                               '4'='Endothelial cells','5'='Fibroblasts','6'='Fibroblasts','7'='Monocytes','8'='Neutrophils',
                               '9'='B cells','10'='DC','11'='NK cells','12'='Fibroblasts','13'='Endothelial cells','14'='Macrophages',
                               '15'='Macrophages','16'='Stromal cells','17'='Macrophages','18'='Endothelial cells','19'='Macrophages',
                               '20'='Fibroblasts','21'='Fibroblasts','22'='B cells')

p1 <- DimPlot(immune.rename1,reduction = "tsne",group.by = "orig.ident")
#DimPlot(immune.combine,reduction = "umap",split.by = "orig.ident")
p2 <- DimPlot(immune.rename1,reduction = "tsne",label = TRUE)
#plot_grid(p1,p2)

immune.rename2 <- RenameIdents(immune.combine,'0'='M_1','1'='F_1','2'='F_2','3'='F_3',
                               '4'='E_1','5'='F_4','6'='F_5','7'='Monocytes','8'='Neutrophils',
                               '9'='B_1','10'='DC','11'='NK','12'='F_6','13'='E_2','14'='M_2',
                               '15'='M_3','16'='S','17'='M_4','18'='E_3','19'='M_5',
                               '20'='F_7','21'='F_8','22'='B_2')
saveRDS("tmp/filter.rename.rds")
