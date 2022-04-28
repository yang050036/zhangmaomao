library(GSVA)
library(limma)
library(msigdbr)
library(BisqueRNA)
library(clusterProfiler)
library(KEGGREST)
library(org.Mm.eg.db)
library(stringr)
data <- readRDS('tmp/filter.rename.rds')
M.data <- subset(data,subset=celltype == 'Macrophages')
Idents(M.data) <- M.data$orig.ident
group <- list(c(c('d3','d5'),'d0'),c(c('d7','d14'),'d0'),c(c('d7','d14'),c('d3','d5')))

markers <- lapply(group,function(x){
  x = unlist(x)
  ident.1 = c(x[1],x[2])
  ident.2 = c(x[3],x[4])
  ident.2 = ident.2[!is.na(ident.2)]
  marker <- FindMarkers(M.data,ident.1=ident.1,ident.2=ident.2,group.by='orig.ident')
})
meta_gene <- read.table('Metabolism_Genes_in_data.xls',header=T,sep='\t')

join <- function(x){
  if(length(x) == 2){
    return(paste(x[1],x[2],sep='&'))
  }else{
    return(x)
  }
}
#i = 1
#group <- list(c(c('d7','d14'),c('d3','d5')))
#group <- list(c(c('d3','d5'),'d0'),c(c('d7','d14'),'d0'))
#group <- list(c(c('d7','d14'),'d0'),c(c('d7','d14'),c('d3','d5')))
Idents(M.data) <- M.data$orig.ident
for(i in 1:length(group)){
  g <- unlist(group[i])
  ident.1 = c(g[1],g[2])
  ident.2 = c(g[3],g[4])
  ident.2 = ident.2[!is.na(ident.2)]
  g1 <- join(ident.1)
  g2 <- join(ident.2)
  title = paste(g1,'vs',g2,sep='_')
  marker <- markers[[i]]
  print(title)
  meta_markers <- marker[intersect(row.names(marker),meta_gene$Gene_name),]
  DE_meta_genes <- row.names(meta_markers)
  xlsx::write.xlsx(meta_markers,file=file.path('10.12/Metabolism/DE',paste(title,".DE.xlsx",sep='')),sheetName='Sheet1')
  
  eg <- bitr(DE_meta_genes, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb='org.Mm.eg.db')
  KEGG <- enrichKEGG(gene = unique(eg$ENTREZID),keyType = "ncbi-geneid",organism = 'mmu',use_internal_data = TRUE,pvalueCutoff = 0.05)
  kegg <- as.data.frame(KEGG)
  #转换ID和name
  kegg$genename <- unlist(lapply(kegg$geneID,  FUN =function(y){
    paste( unlist(lapply(unlist(strsplit(y,"/")),FUN=function(y){eg$SYMBOL[which(eg$ENTREZID ==y)]})) , collapse = "/")} ))
  BP_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="BP",pvalueCutoff = 0.05,readable = TRUE)
  MF_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="MF",pvalueCutoff = 0.05,readable = TRUE)
  CC_go <- enrichGO(gene = unique(eg$ENTREZID),keyType = "ENTREZID",OrgDb = 'org.Mm.eg.db',ont="CC",pvalueCutoff = 0.05,readable = TRUE)
  bp <- as.data.frame(BP_go)
  cc <- as.data.frame(CC_go)
  mf <- as.data.frame(MF_go)
  xlsx::write.xlsx(bp,file=file.path('10.12/Metabolism/GSVA',paste(title,".BP.xlsx",sep='')),sheetName='Sheet1')
  xlsx::write.xlsx(cc,file=file.path('10.12/Metabolism/GSVA',paste(title,".MF.xlsx",sep='')),sheetName='Sheet1')
  xlsx::write.xlsx(mf,file=file.path('10.12/Metabolism/GSVA',paste(title,".CC.xlsx",sep='')),sheetName='Sheet1')
  xlsx::write.xlsx(kegg,file=file.path('10.12/Metabolism/GSVA',paste(title,".KEGG.xlsx",sep='')),sheetName='Sheet1')
  pathlist <- kegg$ID
  kegg_list <- lapply(pathlist,function(x){
    path_info <- keggGet(x)
    genes<-unlist(lapply(path_info[[1]]$GENE,function(x) strsplit(x,';'))) 
    if(!is.null(genes)){
      genelist <- genes[1:length(genes)%%3 ==2]
      return(genelist)
    }else{
      print(x)
    }
  })
  names(kegg_list) <- paste(kegg$ID,kegg$Description,sep=':')
  bpterms <- bp$ID
  ccterms <- cc$ID
  mfterms <- mf$ID
  BP_GO <- get_GO_data("org.Mm.eg.db", "BP", "SYMBOL")
  bp_list <- BP_GO$PATHID2EXTID[bpterms]
  names(bp_list) <- paste(names(BP_GO$PATHID2NAME[bpterms]),BP_GO$PATHID2NAME[bpterms],sep=':')
  print(length(BP_GO$PATHID2EXTID))
  CC_GO <- get_GO_data("org.Mm.eg.db", "CC", "SYMBOL") 
  cc_list <- CC_GO$PATHID2EXTID[ccterms]
  names(cc_list) <-  paste(names(CC_GO$PATHID2NAME[ccterms]),CC_GO$PATHID2NAME[ccterms])
  print(length(CC_GO$PATHID2EXTID))
  MF_GO <- get_GO_data("org.Mm.eg.db", "MF", "SYMBOL")
  mf_list <- MF_GO$PATHID2EXTID[mfterms]
  names(mf_list) <-  paste(names(MF_GO$PATHID2NAME[mfterms]),MF_GO$PATHID2NAME[mfterms])
  print(length(MF_GO$PATHID2EXTID))
  
  if(length(g) == 3){
    group_df <- data.frame(time=g,group = c(g1,g1,g2),row.names = g)
  }else{
    group_df <- data.frame(time=g,group = c(g1,g1,g2,g2),row.names = g)
  }
  sub_data <- SubsetData(M.data,subset.name='orig.ident',accept.value=g)
  sub_data <- sub_data[DE_meta_genes]
  sub_data.ex <- as.SingleCellExperiment(sub_data)
  sub_data.ex$group <- group_df$group[match(sub_data.ex[["orig.ident"]],group_df$time)]
  #exp <- as.data.frame(sub_data@assays$RNA@data)
  #metadata <- sub_data@meta.data
  #exp <- as.matrix(exp)
  sc.eset <- BisqueRNA::SeuratToExpressionSet(sub_data, delimiter="-", position=2, version="v3")
  sc.eset$cellType <- sub_data.ex$celltype
  sc.eset$group <- sub_data.ex$group
  Biobase::exprs(sc.eset) <- as.matrix(sub_data@assays$RNA@data)
  sortdf <- gsva_histgram(sc.eset,kegg_list)
  xlsx::write.xlsx(sortdf,file=file.path('10.12/Metabolism/GSVA',paste(title,".kegg_gsva.xlsx",sep='')),sheetName='Sheet1')
  if(length(sortdf$ID) > 20){sortdf <- sortdf %>% group_by(group) %>% top_n(10,wt=score)}
  h <- gsva_histplot(sortdf)
  height <- (length(sortdf$ID) / 8) * 3.5
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'kegg_11','png',sep='.')),h,height = height)
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'kegg_11','pdf',sep='.')),h,height = height)
  
  sortdf <- gsva_histgram(sc.eset,bp_list)
  xlsx::write.xlsx(sortdf,file=file.path('10.12/Metabolism/GSVA',paste(title,".bp_gsva.xlsx",sep='')),sheetName='Sheet1')
  if(length(sortdf$ID) > 20){sortdf <- sortdf %>% group_by(group) %>% top_n(10,wt=score)}
  h <- gsva_histplot(sortdf)
  height <- (length(sortdf$ID) / 8) * 3.5
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'bp_11','png',sep='.')),h,height = height)
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'bp_11','pdf',sep='.')),h,height = height)
  
  sortdf <- gsva_histgram(sc.eset,cc_list)
  xlsx::write.xlsx(sortdf,file=file.path('10.12/Metabolism/GSVA',paste(title,".cc_gsva.xlsx",sep='')),sheetName='Sheet1')
  if(length(sortdf$ID) > 20){sortdf <- sortdf %>% group_by(group) %>% top_n(10,wt=score)}
  h <- gsva_histplot(sortdf)
  height <- (length(sortdf$ID) / 8) * 3.5
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'cc_11','png',sep='.')),h,height = height)
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'cc_11','pdf',sep='.')),h,height = height)
  
  sortdf <- gsva_histgram(sc.eset,mf_list)
  xlsx::write.xlsx(sortdf,file=file.path('10.12/Metabolism/GSVA',paste(title,".mf_gsva.xlsx",sep='')),sheetName='Sheet1')
  if(length(sortdf$ID) > 20){sortdf <- sortdf %>% group_by(group) %>% top_n(10,wt=score)}
  h <- gsva_histplot(sortdf)
  height <- (length(sortdf$ID) / 8) * 3.5
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'mf_11','png',sep='.')),h,height = height)
  ggsave(file.path('10.12/Metabolism/GSVA',paste(title,'mf_11','pdf',sep='.')),h,height = height)
  i = i + 1
}
