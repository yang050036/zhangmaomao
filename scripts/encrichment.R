library(Seurat)
library(ggplot2)
library(KEGG.db)
##GO/KEGG
data<-readRDS("/mnt/beegfs/Research/Users/yangyupeng/scRNA/zhangmaomao/seurat/tmp/filter.rename.rds")##文件
data$cluster <- Idents(data)##加一列（属于哪个cluster）
data1<-subset(data,subset = celltype=="B cells")##提取巨噬细胞
marker01<- FindMarkers(data1, ident.1 ="d0",ident.2 = "d1",group.by = "orig.ident")
marker<- bitr(rownames(marker01), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
GO_BP<- enrichGO(gene = marker$ENTREZID, OrgDb="org.Mm.eg.db", keyType = "ENTREZID",ont = "BP", readable= TRUE)
GO_CC<- enrichGO(gene = marker$ENTREZID, OrgDb="org.Mm.eg.db", keyType = "ENTREZID",ont = "CC", readable= TRUE)
GO_MF<- enrichGO(gene = marker$ENTREZID, OrgDb="org.Mm.eg.db", keyType = "ENTREZID",ont = "MF", readable= TRUE)
B<-barplot(GO_BP, showCategory=30,title="EnrichmentGO_BP day0 vs day1") 
C<-barplot(GO_CC, showCategory=30,title="EnrichmentGO_CC day0 vs day1")
M<-barplot(GO_MF, showCategory=30,title="EnrichmentGO_MF day0 vs day1")

marker01<- FindMarkers(data1, ident.1 ="d0",ident.2 = "d1",group.by = "orig.ident")
marker<- bitr(rownames(marker01), fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Mm.eg.db")
KEGG<-enrichKEGG(gene = marker$ENTREZID, organism = "mmu", keyType = "kegg",use_internal_data = TRUE)
K<-barplot(KEGG, showCategory=30,title="EnrichmentKEGG day0 vs day1")
