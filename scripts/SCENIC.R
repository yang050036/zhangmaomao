#!/usr/bin/Rscript

library(optparse)
optionlist <- list(
  make_option('--input',action='store_true',type='character',default=FALSE,
              help="Please input a database from rds."),
  make_option('--output',action='store_true',type='character',default=FALSE,
              help="Output path."),
  make_option('--gene_list',action='store_true',type='character',default=FALSE,
              help="gene list file."),
  make_option('--cell_type',action='store_true',type='character',default=FALSE,
              help="cell type."),
  make_option('--resource',action='store_true',type='character',default=FALSE,
              help="resource path way.")
)  
opt <- parse_args(OptionParser(option_list=optionlist))

library(SCENIC)
library(RcisTarget)
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
library(gtools)
library(parallel)
library(doParallel)
library(loomR)
library(ComplexHeatmap)

if(dir.exists(opt$output)){
  setwd(opt$output)
}else{
  dir.create(opt$output)
  setwd(opt$output)
}

if(!dir.exists('int')){
  dir.create('int')
}

#day <- readRDS('tmp/merge.rds')
rds <- readRDS(opt$input)
sub_rds <- subset(rds, subset = celltype == opt$cell_type)
DefaultAssay(sub_rds) <- "RNA"
exprMat<-GetAssayData(object = sub_rds,slot = "counts")
cellInfo <- sub_rds@meta.data
colnames(cellInfo) <- c('time','nGene','nUMI','percent.mt','integreted','seurat_clusters','cellType1','filter','celltype','cluster')
cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "Class"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "CellType"
saveRDS(cellInfo, file="int/cellInfo.Rds")
org <- "mgi" 
dbDir <- opt$resource
myDatasetTitle <- "musscRNA" # Create a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20)
saveRDS(scenicOptions, file="scenicOptions.Rds")

#filter genes without of database
markers <- read.table(opt$gene_list,header=T,sep=' ')

top20 <- markers %>% group_by(cluster) %>% top_n(n=100,wt=avg_logFC)
genelist <- unique(markers$gene)
dbFilePath <- getDatabases(scenicOptions)[[1]]
motifRankings <- importRankings(dbFilePath)
genesInDatabase <- colnames(getRanking(motifRankings))
genesLeft_minCells<-rownames(exprMat)
length(genesLeft_minCells)#19025
genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)#14538
genesKept <- genesLeft_minCells_inDatabases
exprMat_filtered <- exprMat[genelist, ]

#共表达网络
runCorrelation(as.matrix(exprMat_filtered), scenicOptions)
exprMat_filtered_log <- log2(as.matrix(exprMat_filtered)+1)
runGenie3(as.matrix(exprMat_filtered_log), scenicOptions)
save(exprMat_filtered,scenicOptions,file = paste(opt$cell_type,"_GENIE3_data.Rdata",sep=''))

#load(paste(opt$cell_type,"_GENIE3_data.Rdata",sep=''))
cellInfo <- ('int/cellInfo.Rds')
#exprMat<-exprMat_filtered
exprMat_log <- log2(as.matrix(exprMat_filtered)+1)
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 50
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
scenicOptions@settings$nCores <- 1
runSCENIC_3_scoreCells(scenicOptions, as.matrix(exprMat_log))
scenicOptions@fileNames$output["loomFile",] <- "output/auc.loom"
fileName <- getOutName(scenicOptions, "loomFile")
if (file.exists(fileName)){
  file.remove(fileName)
}
export2loom <- export2scope <- function(scenicOptions, dgem, hierarchy=c("SCENIC", "", ""), 
                                        addAllTsnes=TRUE)
{
  regulonType="Motif" # TODO
  # TODO: what about if there is no dgem, but only normalized?
  
  if(length(hierarchy) > 3) stop("Maximum three hierarchy levels")
  if(length(hierarchy) < 3) hierarchy <- c(hierarchy, rep("", 5-length(hierarchy)))
  fileName <- getOutName(scenicOptions, "loomFile")
  if(file.exists(fileName)) stop("File '", fileName, "' already exists.")
  
  # TODO: ask max about order of samples tsne-expr-info
  suppressPackageStartupMessages(require(SCopeLoomR))
  
  # Default embedding (e.g. t-SNE or PCA coordinates)
  if(!file.exists(tsneFileName(scenicOptions))) stop(paste("Default 2D projection is not available:", tsneFileName(scenicOptions)))
  defaultTsne <- readRDS(tsneFileName(scenicOptions))
  defaultTsne_name <- paste("SCENIC t-SNE:", defaultTsne$type)
  defaultTsne <- defaultTsne$Y
  cellOrder <- rownames(defaultTsne)
  if(!all(colnames(dgem) == cellOrder)) warning("tSNE and matrix cell names do not match.")
  
  # Cell info:
  cellInfo <- loadFile(scenicOptions, getDatasetInfo(scenicOptions, "cellInfo"), ifNotExists="null")
  if(!all(rownames(cellInfo) == cellOrder)) warning("CellInfo cell names (order?) do not match.")
  for(cn in colnames(cellInfo)) if(is.numeric(cellInfo[,cn]) & is.character(cellInfo[,cn])) stop(paste0(cn, "should be either numeric or character/factor."))
  if(getSettings(scenicOptions, "verbose")) {
    message("The folowing cell metadata will be added:")
    print(data.frame(type=cbind(sapply(cellInfo, class))))
  }
  
  # Start adding...
  fileName <- getOutName(scenicOptions, "loomFile")
  loom <- build_loom(file.name=fileName,
                     dgem=dgem,
                     title=getDatasetInfo(scenicOptions,"datasetTitle"),
                     genome=getDatasetInfo(scenicOptions,"org"), # Just for user information, not used internally
                     default.embedding=defaultTsne,
                     default.embedding.name=defaultTsne_name)
  
  add_hierarchy(loom=loom, hierarchy=create_hierarchy(level.1.name=hierarchy[1], 
                                                      level.2.name=hierarchy[2], 
                                                      level.3.name=hierarchy[3]))
  
  # Known cell information/annotation  
  if(!is.null(cellInfo))
  {
    cellInfo <- data.frame(cellInfo)
    cellInfo <- cellInfo[,colnames(cellInfo) != "nGene", drop=FALSE]
    cellInfo <- cellInfo[,colnames(cellInfo) != "nUMI", drop=FALSE]
    
    # Add annotation
    for(cn in colnames(cellInfo))
    {
      isMetric <- FALSE
      isAnnotation <- FALSE
      if(is.character(cellInfo[,cn]) || is.factor(cellInfo[,cn])) 
      {
        isAnnotation <- TRUE
        cellInfo[,cn] <- as.character(cellInfo[,cn])
      }else{
        if(all(!is.na(as.numeric(as.character(cellInfo[,cn]))))) isMetric <- TRUE
      }
      add_col_attr(loom=loom, key=cn, value=cellInfo[,cn], as.annotation=isAnnotation, as.metric=isMetric)
    }
  }
  
  # Regulons AUC matrix
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  if(!all(colnames(regulonAUC) == cellOrder)) warning("regulonAUC cell names (order?) do not match.")
  add_scenic_regulons_auc_matrix(loom=loom, regulons.AUC=AUCell::getAUC(regulonAUC))
  
  # Regulons (gene list)
  regulons <- loadInt(scenicOptions, "aucell_regulons")
  motifEnrichment <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
  motifEnrichment <- motifEnrichment[grep("transfac_pro__", motifEnrichment[["motif"]], invert = T),] #TODO
  regulonThresholds <- loadInt(scenicOptions, "aucell_thresholds")
  if(!"aucThr" %in% names(regulonThresholds[[1]]))
  {
    # regulonThresholds <- lapply(regulonThresholds, function(x) list(selected=x))
    warning("The binarized regulon activity will not been added to the loom.")
    regulonThresholds <- NULL
  }
  add_scenic_regulons(loom=loom
                      , regulons=regulons
                      , regulon.threshold.assignments=regulonThresholds # Optional
                      , regulon.enrichment.table=motifEnrichment # Optional
                      
  )
  
  # # Alternative t-SNE #TODO
  # load(paste0(scenicDir, "int/3.3_tsneRegulonAUC_PCA_50pc50perpl.RData"))
  # add_embedding(loom=loom, embedding=tSNE$Y, name="SCENIC (t-SNE on AUC)")
  if(addAllTsnes)
  {
    tsnePath <- dirname(getSettings(scenicOptions, "tSNE_filePrefix"))
    allTsneFiles <- grep(pattern = paste0(basename(getSettings(scenicOptions, "tSNE_filePrefix")), ".*\\.Rds")
                         , x=list.files(tsnePath)
                         , value = T)
    
    for(tsneFile in allTsneFiles)
    {
      tryCatch(
        add_embedding(loom=loom, embedding=readRDS(file.path(tsnePath, tsneFile))$Y, name=paste0("SCENIC: ",  gsub(".Rds","", tsneFile, fixed=TRUE)))
        , error=function(e) print(paste0("Cannot add tsne: ", e$message)))
    }
  }
  
  finalize(loom=loom)
  if(getSettings(scenicOptions, "verbose")) message("Loom file saved as:\t", fileName)
}
#export2loom(scenicOptions, exprMat)
#aucellApp <- plotTsne_AUCellApp(scenicOptions, as.matrix(exprMat_log))
#savedSelections <- shiny::runApp(aucellApp)
#runSCENIC_4_aucell_binarize(scenicOptions)

#tSNE结果
#nPcs <- c(5)
#scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
#scenicOptions@settings$seed <- 50
#fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
#fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
#fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
#fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
#par(mfrow=c(length(nPcs), 4))
#plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="celltype", cex=.5)
#scenicOptions@settings$defaultTsne$aucType <- "AUC"
#scenicOptions@settings$defaultTsne$dims <- 5
#scenicOptions@settings$defaultTsne$perpl <- 15
#saveRDS(scenicOptions, file="int/scenicOptions.Rds")

#motifsAUC <- loadInt(scenicOptions, "motifs_AUC")
#motifsAUC <- motifsAUC[row.names(motifsAUC),]
#motifsActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$celltype),
#                                     function(cells) rowMeans(getAUC(motifsAUC)[,cells]))

#regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

#h <- pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3, 
#                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
#                   treeheight_row=10, treeheight_col=10, border_color=NA)
#ggsave('output/known_cell_Regulators.pdf',h)
#ggsave('output/known_cell_Regulators.png',h)
#cell_type=opt$cell_type
#sub_rds <- subset(rds, subset = celltype == cell_type)
#sub_cluster <- mixedsort(unique(sub_rds$cluster))
#cellInfo <- subset(cellInfo,subset=cluster %in% M)
motifsAUC <- loadInt(scenicOptions, "motifs_AUC")
motifEnrichmentTable <- addMotifAnnotation(motifsAUC, nesThreshold=3,
                                           motifAnnot=motifAnnotations_mgi)
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes")
anotatedTfs  <- lapply(split(motifEnrichment_selfMotifs_wGenes$TF_highConf,
                             motifEnrichmentTable$geneSet),
                       function(x) {
                         genes <- gsub(" \\(.*\\). ", "; ", x, fixed=FALSE)
                         genesSplit <- unique(unlist(strsplit(genes, "; ")))
                         return(genesSplit)
                       })

data <- readRDS('../tmp/filter.rename.rds')
celltype_list <- unique(data$celltype)
work_path <- getwd()
lapply(celltype_list,function(x){
  x <- gsub(' ','_',x)
  setwd(file.path(work_path,x))
  scenicOptions <- readRDS(file.path('int','scenicOptions.Rds'))
  cellInfo<- readRDS(file.path('int','cellInfo.Rds'))
  cellInfo$cluster <- factor(cellInfo$cluster, levels = mixedsort(unique(cellInfo$cluster)))
  regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
  regulonActivity_byCellType <- getAUC(regulonAUC)
  regulonActivity_byCellType <- regulonActivity_byCellType[!grepl('extended',row.names(regulonActivity_byCellType)),]
  regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
  
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]}
  color <- gg_color_hue(length(mixedsort(unique(cellInfo$cluster))))
  names(color) <- mixedsort(unique(cellInfo$cluster))
  
  time_color <- c("d0" = "#020101FF","d1" = "#BF747EFF","d3" = "#415F49FF","d5" = "#3A173EFF","d7" = "#7D908FFF","d14" = "#C2E95EFF")
  time_keep <- mixedsort(unique(cellInfo$time))
  time_color <- time_color[time_keep]
  cluster_order <- cellInfo$time[mixedorder(cellInfo$cluster)]
  time_order <- lapply(unique(mixedsort(cellInfo$cluster)),function(cluster){
    list <- which((mixedsort(cellInfo$cluster) == cluster)==TRUE)
    order <- list[mixedorder(cluster_order[list])]
    return(order)
  })
  colorder <- mixedorder(cellInfo$cluster)[unlist(time_order)]
  cluster_anno <- HeatmapAnnotation(cluster = factor(cellInfo$cluster[colorder],levels=mixedsort(unique(cellInfo$cluster))),
                                    time = factor(cellInfo$time[colorder],levels=mixedsort(unique(cellInfo$time))),
                                    col=list(cluster=color,time=time_color),show_legend = c("cluster" = TRUE,"time" = TRUE))
  
  #timeorder <- mixedorder(cellInfo$time)
  mtx <- regulonActivity_byCellType_Scaled

  #mtx <- mtx[,timeorder]
  mtx <- mtx[,colorder]
  heatmap_legend_param <- list(grid_height = unit(8, "mm"))
  print(file.path('output',paste(x,'heatmap','png',sep='.')))
  png(file.path('output',paste(x,'heatmap','png',sep='.')),width=900,height = 1000,res=150)
  print(Heatmap(mtx,top_annotation = cluster_anno,cluster_rows = T,cluster_columns = F,show_column_names = F,column_title = NULL,
          row_names_gp = gpar(fontsize = 7),heatmap_legend_param=heatmap_legend_param,column_split = mixedsort(cellInfo$cluster)))
  dev.off()
  pdf(file.path('output',paste(x,'heatmap','pdf',sep='.')))
  print(Heatmap(mtx,top_annotation = cluster_anno,cluster_rows = T,cluster_columns = F,show_column_names = F,column_title = NULL,
          row_names_gp = gpar(fontsize = 7),heatmap_legend_param=heatmap_legend_param,column_split = mixedsort(cellInfo$cluster)))
  dev.off()
})
#auc <- as.data.frame(scale(t(regulonActivity_byCellType), center = T, scale=T))
#auc$cluster <- cellInfo$cluster
#topRegulators <- reshape2::melt(auc)

#h <- pheatmap(motifsActivity_bycluster_Scaled, #fontsize_row=3, 
#              color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-3, 3, length.out = 100),
#              treeheight_row=10, treeheight_col=10, border_color=NA)
#ggsave(file.path('output',paste(cell_type,'heatmap','pdf',sep = '.')),h)
#ggsave(file.path('output',paste(cell_type,'heatmap','png',sep = '.')),h)

#导入seurat的结果中
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
#auc_loom <- as.loom(AUCmatrix, filename = 'int/auc.loom',verbose=T)
#auc_loom$close_all()
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(rds, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'tmp/scRNAauc.rds')

##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(rds, BINmatrix)
scRNAbin@assays$integrated <- NULL
saveRDS(scRNAbin, 'tmp/scRNAbin.rds')

lapply(RegulonName_BIN,function(x){
  p <- FeaturePlot(scRNAbin,features=x,reduction = 'tsne',slot='scale.data',blend = FALSE)
  ggsave(file.path('picture/scenic_seurat',paste(x,'.point.png')),p,width = 10)
  ggsave(file.path('picture/scenic_seurat',paste(x,'.point.pdf')),p,width = 10)
})

#newwork
signifMotifNames <- motifEnrichmentTable$motif[1:3]

incidenceMatrix <- getSignificantGenes(geneLists$hypoxia, 
                                       motifRankings,
                                       signifRankingNames=signifMotifNames,
                                       plotCurve=TRUE, maxRank=5000, 
                                       genesFormat="incidMatrix",
                                       method="aprox")$incidMatrix

library(reshape2)
.libPaths('/mnt/beegfs/Research/Software/10X/R-4.0.0/lib64/R/library')
edges <- melt(incidenceMatrix)
edges <- edges[which(edges[,3]==1),1:2]
colnames(edges) <- c("from","to")
library(visNetwork)
motifs <- unique(as.character(edges[,1]))
genes <- unique(as.character(edges[,2]))
nodes <- data.frame(id=c(motifs, genes),   
                    label=c(motifs, genes),    
                    title=c(motifs, genes), # tooltip 
                    shape=c(rep("diamond", length(motifs)), rep("elypse", length(genes))),
                    color=c(rep("purple", length(motifs)), rep("skyblue", length(genes))))
visNetwork(nodes, edges) %>% visOptions(highlightNearest = TRUE, 
                                        nodesIdSelection = TRUE)

#rss
library(scFunctions)
library(tidyverse)
library(ggridges)
number_of_regulon_clusters <- 10
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
kmeans_thresholds <- auc_thresh_kmeans(regulonAUC)
binary_regulons <- binarize_regulons(regulonAUC,kmeans_thresholds)
joined_bin_reg <- binary_regulons %>% reduce(left_join,by="cells")
row.names(joined_bin_reg) <- joined_bin_reg$cells
joined_bin_reg <- joined_bin_reg[2:ncol(joined_bin_reg)]
binary_regulons_trans <- as.matrix(t(joined_bin_reg))
metadata_sub <- readRDS('int/cellInfo.Rds')
metadata_sub$cluster <- factor(metadata_sub$cluster,levels=mixedsort(unique(metadata_sub$cluster)))
rrs_df <- calculate_rrs(metadata_sub, binary_regulons = binary_regulons_trans,cell_type_column="cluster")

calculate_csi <- function (regulonAUC, calc_extended = FALSE, verbose = FALSE) 
{
  compare_pcc <- function(vector_of_pcc, pcc) {
    pcc_larger <- length(vector_of_pcc[vector_of_pcc > pcc])
    if (pcc_larger == length(vector_of_pcc)) {
      return(0)
    }
    else {
      return(length(vector_of_pcc))
    }
  }
  calc_csi <- function(reg, reg2, pearson_cor) {
    test_cor <- pearson_cor[reg, reg2]
    total_n <- ncol(pearson_cor)
    pearson_cor_sub <- subset(pearson_cor, rownames(pearson_cor) == 
                                reg | rownames(pearson_cor) == reg2)
    sums <- apply(pearson_cor_sub, MARGIN = 2, FUN = compare_pcc, 
                  pcc = test_cor)
    fraction_lower <- length(sums[sums == nrow(pearson_cor_sub)])/total_n
    return(fraction_lower)
  }
  regulonAUC_sub <- regulonAUC@assays@data@listData$AUC
  if (calc_extended == TRUE) {
    regulonAUC_sub <- subset(regulonAUC_sub, grepl("extended", 
                                                   rownames(regulonAUC_sub)))
  }
  else if (calc_extended == FALSE) {
    regulonAUC_sub <- subset(regulonAUC_sub, !grepl("extended", 
                                                    rownames(regulonAUC_sub)))
  }
  regulonAUC_sub <- t(regulonAUC_sub)
  pearson_cor <- cor(regulonAUC_sub)
  pearson_cor_df <- as.data.frame(pearson_cor)
  pearson_cor_df$regulon_1 <- rownames(pearson_cor_df)
  pearson_cor_long <- pearson_cor_df %>% gather(regulon_2, 
                                                pcc, -regulon_1) %>% mutate(regulon_pair = paste(regulon_1, 
                                                                                                 regulon_2, sep = "_"))
  regulon_names <- unique(colnames(pearson_cor))
  num_of_calculations <- length(regulon_names) * length(regulon_names)
  csi_regulons <- data.frame(matrix(nrow = num_of_calculations, 
                                    ncol = 3))
  colnames(csi_regulons) <- c("regulon_1", "regulon_2", "CSI")
  num_regulons <- length(regulon_names)
  f <- 0
  for (reg in regulon_names) {
    if (verbose == TRUE) {
      print(reg)
    }
    for (reg2 in regulon_names) {
      f <- f + 1
      fraction_lower <- calc_csi(reg, reg2, pearson_cor)
      csi_regulons[f, ] <- c(reg, reg2, fraction_lower)
    }
  }
  csi_regulons$CSI <- as.numeric(csi_regulons$CSI)
  return(csi_regulons)
}
regulons_csi <- calculate_csi(regulonAUC,  calc_extended = FALSE)

####################
rrs_df <- read.table('NK.rss.tsv',sep='\t',header = T)
rrs_df_wide <- rrs_df %>% spread(cell_type,RSS)
rownames(rrs_df_wide) <- rrs_df_wide$regulon
#rrs_df_wide <- rrs_df_wide[,2:ncol(rrs_df_wide)]
rrs_df_wide <- data.frame('NK' = rrs_df_wide[,2:ncol(rrs_df_wide)],row.names = row.names(rrs_df_wide))
h <- pheatmap(rrs_df_wide,border_color=NA,
         show_colnames = TRUE,
         color = viridis(n = 10),
         cluster_cols = F,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         fontsize_row = 5,cellwidth = 20)
png('NK.rss.png',width = 3000,height = 3000,res=300)
h
dev.off()
pdf('NK.rss.pdf')
h
dev.off()
#rrs_df_wide_specific <- rrs_df_wide[apply(rrs_df_wide,MARGIN = 1 ,FUN =  function(x) any(x > 0.4)),]
#rrs_df_nona <- subset(rrs_df,RSS > 0)
rank_plot <-plot_rrs_ranking(rrs_df,
                 "RMPH",
                 ggrepel_force = 1,
                 ggrepel_point_padding = 0.2,
                 top_genes = 4,
                 plot_extended = FALSE)

png('NK.rank.png',width = 3000,height = 3000,res=300)
rank_plot
dev.off()
pdf('NK.rank.pdf')
rank_plot
dev.off()
#csi

regulons_csi <- read.table('NK.csi.tsv',sep='\t',header = T)

h <- plot_csi_modules(regulons_csi,
                 nclust = 10,
                 font_size_regulons = 8)

png('NK.csi.png',width = 3000,height = 3000,res=300)
h
dev.off()
pdf('NK.csi.pdf')
h
dev.off()


