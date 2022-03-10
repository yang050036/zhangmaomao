library(Seurat)
library(ggplot2)
library(magrittr)
library(dplyr)
library(cowplot)
library(metap)
library(egg)
library(reshape2)
library(RColorBrewer)
library(SingleCellExperiment)
library(pheatmap)

#read file from cellranger of 10X genomeics
files <- list.files(opt$files_path,all.files = TRUE,no..=TRUE)
filelist <- lapply(files,function(file){
  project=file
  file = file.path(opt$files_path,file)
  data = Read10X(file)
  if(is.list(data)){
    data = data$`Gene Expression`
  }
  CreateSeuratObject(counts = data,
                     project = project,min.cells = 3,min.features = 200)
})
files_num <- length(files)
file_other <- c()
for (i in 2:files_num){
  file_other[[i-1]] <- filelist[[i]]
}
#merge data
data <- merge(filelist[[1]],y = file_other,add.cell.ids = files)
data[['percent.mt']] <- PercentageFeatureSet(data,pattern = "^MT")
data <- subset(data,subset = nFeature_RNA >200 & nFeature_RNA <2500 & percent.mt<30)

#normalize data for each sample and integrated data with CCA method
data.list <- SplitObject(data,split.by = 'orig.ident')
data.list <- lapply(X = data.list,FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x,selection.method = 'vst',nfeatures = 2000)
})
immune.anchors <- FindIntegrationAnchors(object.list = data.list,dims = 1:20)
data <- IntegrateData(anchorset = immune.anchors,dims = 1:20)
DefaultAssay(data) <- "integrated"

#scale data with all variable features, dimension reduction and cluster cells. 
all.genes <- row.names(data)
data <- ScaleData(data,features = all.genes,verbose = FALSE)
data <- RunPCA(data,npcs = 30,verbose = FALSE)
data <- RunTSNE(data,reduction = "pca",dims = 1:20)
data <- FindNeighbors(data,reduction = 'pca',dims = 1:20)
data <- FindClusters(data,resolution = 0.5)

#defined cluster with celltype
data = RenameIdents(data,
                    "0" = "Neutrophils",
                    "1" = "Neutrophils",
                    "2" = "Monocytes",
                    "3" = "Monocytes",
                    "4" = "Neutrophils",
                    "5" = "NK cells",
                    "6" = "Neutrophils",
                    "7" = "Monocytes",
                    "8" = "B cells",
                    "9" = "Monocytes",
                    "10" = "Neutrophils",
                    "11" = "Neutrophils",
                    "12" = "Neutrophils",
                    "13" = "T cells",
                    "14" = "B cells",
                    "15" = "Neutrophils",
                    "16" = "Monocytes",
                    "17" = "Monocytes")
data$celltype = Idents(data)
saveRDS(data, "celltype.rds")
