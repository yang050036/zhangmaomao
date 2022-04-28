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
library(htmlwidgets)
library(rmarkdown)
library(gtools)


seurat_rds <- readRDS('tmp/filter.rename.rds')
DefaultAssay(seurat_rds) <- "RNA"
sub_rds <- subset(seurat_rds,subset=celltype=='Macrophages')
sub_rds$times <- as.numeric(unlist(lapply(sub_rds$orig.ident,function(x){gsub('d','',x)})))
sub_rds$orig.ident <- factor(sub_rds$orig.ident,levels = mixedsort(unique(sub_rds$orig.ident)))
sub_rds$rename <- factor(sub_rds$cluster,levels = mixedsort(unique(sub_rds$cluster)))
data <- sub_rds@assays$RNA@data
pd <- sub_rds@meta.data
new_pd<-select(pd,orig.ident,nCount_RNA,nFeature_RNA,percent.mt,times,rename)
#new_pd$celltype <- Idents(sub_rds)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
cds <- new_cell_data_set(data,cell_metadata  = new_pd,gene_metadata  = fData)
cds <- preprocess_cds(cds, num_dim = 30,norm_method = "log")
cds = reduce_dimension(cds, reduction_method="UMAP")
cds <- cluster_cells(cds,restion = 0.5,reduction_method = "UMAP")
#p <- plot_cells(cds,color_cells_by = "cluster",reduction_method = "UMAP")

cds <- learn_graph(cds)
get_earliest_principal_node <- function(cds, group.by,time_bin="d0"){
  cell_ids <- which(colData(cds)[, group.by] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
time_bin <- get_earliest_principal_node(cds,group.by = 'times',time_bin = '0')
sample_bin <- get_earliest_principal_node(cds,group.by = 'orig.ident',time_bin = 'd0')
cluster_bin <- get_earliest_principal_node(cds,group.by = 'rename',time_bin = 'M_1')
time <- order_cells(cds, root_pr_nodes=time_bin)
sample <- order_cells(cds, root_pr_nodes=sample_bin)
cluster <- order_cells(cds, root_pr_nodes=cluster_bin)
time_ti <- plot_cells(time,
                      color_cells_by = "times",
                      label_cell_groups=FALSE,
                      label_leaves=FALSE,
                      label_branch_points=FALSE,
                      graph_label_size=1.5)
sample_ti <- plot_cells(sample,
                        color_cells_by = "orig.ident",
                        label_cell_groups=FALSE,
                        label_leaves=FALSE,
                        label_branch_points=FALSE,
                        graph_label_size=1.5)
cluster_ti <- plot_cells(cluster,
                         color_cells_by = "rename",
                         label_cell_groups=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=FALSE,
                         graph_label_size=1.5)
#dir.create('picture/ti')
ggsave(file.path('picture/ti','time_ti.png'),time_ti)
ggsave(file.path('picture/ti','time_ti.svg'),time_ti)
ggsave(file.path('picture/ti','sample_ti.png'),sample_ti)
ggsave(file.path('picture/ti','sample_ti.svg'),sample_ti)
ggsave(file.path('picture/ti','cluster_ti.png'),cluster_ti)
ggsave(file.path('picture/ti','cluster_ti.svg'),cluster_ti)
#3D图
#cds_3d <- reduce_dimension(cds, max_components = 3)
#cds_3d <- cluster_cells(cds_3d)
#cds_3d <- learn_graph(cds_3d)
#cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
#cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="partition")
#saveWidget(cds_3d_plot_obj,'/mnt/beegfs/Research/Users/yangyupeng/scRNA/zhangmaomao/seurat/picture/ti3.html')

#基因层面探索拟时序分析
diff_gene <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
id <- row.names(subset(diff_gene, q_value < 0.05))
CASE_genes <- c("Xkr4", "Sox17", "Lypla1",'Tectb')
CASE_lineage_cds <- cds[rowData(cds)$gene_short_name %in% CASE_genes,colData(cds)$celltype %in% c("M_1")]
ti2 <- plot_genes_in_pseudotime(CASE_lineage_cds,color_cells_by="celltype")

pdf('picture/ti_gene')
ti2
dev.off()

#ggsave('picture/ti_1.png',ti1)
#ggsave('picture/ti_gene.png',ti2)

##################SCORPIUS软件###########
library(SCORPIUS)
rds <- readRDS('tmp/filter.rename.rds')
M_rds <- subset(rds,subset=celltype=='Macrophages')
M_rds$orig.ident <- factor(M_rds$orig.ident,levels = mixedsort(as.character(unique(M_rds$orig.ident))))
M_rds$celltype <- factor(M_rds$celltype,levels = mixedsort(unique(M_rds$celltype)))
M_rds$cluster <- factor(M_rds$cluster,levels = mixedsort(as.character(unique(M_rds$cluster))))
M_rds$times <- as.numeric(unlist(lapply(M_rds$orig.ident,function(x){gsub('d','',x)})))
#DefaultAssay(M_rds) <- 'integrated'
#expressed_genes <- row.names(M_rds)
top20 <- M_markes %>% group_by(cluster) %>% top_n(n=20,wt = avg_logFC)

mex <- M_rds@assays$RNA@data[genes,]
expression <- t(as.matrix(mex))
groupname <- M_rds$orig.ident
clustername <- Idents(M_rds)
space <- reduce_dimensionality(expression,dist='spearman',ndim=30)#spearman,pearson,euclidean,cosine,manhattan
traj <- infer_trajectory(space,k=9)
p <- mydraw_trajectory_plot(space,progression_group=groupname,path=traj$path,contour=FALSE,point_size=0.5)


ggsave('picture/ti/SCORPIUS_trajectory.png',p)
ggsave('picture/ti/SCORPIUS_trajectory.svg',p)

#对draw_trajectory_heatmap函数进行了edit，将annotation的名字改成了cluster
gimp <- gene_importances(expression,traj$time,num_permutations = 0, num_threads = 8)
gene_sel <- gimp[1:50,]
expr_sel <- expression[,gene_sel$gene]
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
h <- newheatmap(expr_sel, traj$time, clustername, modules,show_labels_row=TRUE,annotation='Sample')

ggsave('picture/ti/SCORPIUS_heatmap.png',h)
ggsave('picture/ti/SCORPIUS_heatmap.svg',h)
xlsx::write.xlsx(as.data.frame(modules),file=file.path('picture/ti',paste("SCORPIUS",".modules.xlsx",sep='')),sheetName='Sheet1',row.names=F)

top4 <- modules %>% group_by(module) %>% top_n(n=4,wt = within_module_ordering)
plots <- list()
for(i in top4$feature){
  p <- draw_trajectory_plot(space,progression_group=expression[,which(colnames(expression) == i)],progression_group_palette=c("lightgrey", "blue"),point_size=0.5) +
    ggtitle(i) + theme(legend.title = element_blank())
  plots[[i]] = p
}

ggsave("picture/ti/SCORPIUS_trajectory_gene.png",ggarrange(plots=plots,ncol=4),width=30,height=40,limitsize=FALSE)
ggsave("picture/ti/SCORPIUS_trajectory_gene.pdf",ggarrange(plots=plots,ncol=4),width=30,height=40,limitsize=FALSE,onefile=F)

