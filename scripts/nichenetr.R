#########受配体#####
library(nichenetr)
library(Seurat)
library(tidyverse)
data <- readRDS('tmp/rename.final.rds')
ligand_target_matrix = readRDS("ligand_target_matrix.rds")
lr_network = readRDS("lr_network.rds")
weighted_networks = readRDS("weighted_networks.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
## receiver
receiver = "Fibroblasts"
Idents(data) <- data$celltype
expressed_genes_receiver = get_expressed_genes(receiver, data, pct = 0.10)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
## sender
sender_celltypes = c("Macrophages")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, data, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

seurat_obj_receiver= subset(data, idents = receiver)
seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[["celltype"]])
DE_table_receiver = FindMarkers(object = data, ident.1 = "Macrophages", ident.2 = "Fibroblasts", min.pct = 0.10) %>% rownames_to_column("gene")

geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) >= 0.25) %>% pull(gene)
#geneset_oi = row.names(seurat_obj_receiver) %>% .[. %in% rownames(ligand_target_matrix)]
geneset_oi =geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(dplyr::desc(pearson)))

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

xlsx::write.xlsx(as.data.frame(ligand_activities),file=file.path('10.12/la/4',"Macrophages_to_Fibroblasts.ligand_activities.xlsx"),sheetName='Sheet1',row.names=F)


#以下数据均取前top20的配体
#活性靶基因预测
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links)
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands from sender cell\n(Macrophages)","Predicted target genes in receiver cells(Fibroblasts)", color = "purple",legend_position = "right", 
                      x_axis_position = "top",legend_title = "Regulatory potential")  + 
  theme(axis.text.x = element_text(face = "italic")) + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
ggsave(file.path('10.12/la/4',paste('Macrophages_to_Fibroblasts','Active_target_gene_inference','png',sep='.')),p_ligand_target_network,width=15)
ggsave(file.path('10.12/la/4',paste('Macrophages_to_Fibroblasts','Active_target_gene_inference','pdf',sep='.')),p_ligand_target_network,width=15)

#排名靠前的配体的受体，只考虑了文献和公开数据库中记载的配体-受体相互作用之后
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()
p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% 
  t() %>% 
  make_heatmap_ggplot("Ligands from sender cells\n(Macrophages)","Receptors in receiver cells\n(Fibroblasts)", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)",legend_position = "right")
