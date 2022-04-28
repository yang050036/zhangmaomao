#!/usr/bin/Rscript

library(optparse)
optionlist <- list(
  make_option('--file',action='store_true',type='character',default=FALSE,
              help="Please input a database from rds."),
  make_option('--project',action='store_true',type='character',default=FALSE,
              help="Providing a name of project."),
  make_option('--output',action='store_true',type='character',default=FALSE,
              help="Output path."),
  make_option('--species',action='store_true',type='character',default=FALSE,
              help="Species of you data"),
  make_option('--font_path',action='store_true',type='character',default=FALSE,
              help="There need a font for cellchat sortware.")
)  
opt <- parse_args(OptionParser(option_list=optionlist))

library(CellChat)
library(ggplot2)
library(Seurat)
library(showtext)
library(ggalluvial)
options(stringsAsFactors = FALSE)

tmp = file.path(opt$output,'tmp')
communication = file.path(opt$output,'communication')
project_path = file.path(communication,opt$project)
file_list = c(opt$output,communication,project_path)
for(file in file_list){
  if(!dir.exists(file)){
    dir.create(file)
  }
}
font_add("Arial", opt$font_path)
#Arial <- Type1Font(family = "Arial",metrics = c("/mnt/beegfs/Research/Pipeline/10X/pipeline/scripts/bin/R/Arial.afm"))
#pdfFonts(Arial = Arial)
rds <- readRDS(opt$file)
data.input <- rds@assays$RNA@data
identity <- data.frame(group=rds$celltype,row.names=names(rds$celltype))
cellchat <- createCellChat(data=data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
if(opt$species == 'mouse'){
  CellChatDB <- CellChatDB.mouse
}else if(opt$species == 'human'){
  CellChatDB <- CellChatDB.human
}
#CellchatDB.use <- subsetDB(CellChatDB,search='Secreted Signaling')#侧重方向包含3个"Secreted Signaling","ECM-Receptor","Cell-Cell Contact"
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers = 4)#多线程
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
if(opt$species == 'mouse'){
  cellchat <- projectData(cellchat, PPI.mouse)
}else if(opt$species == 'human'){
  cellchat <- projectData(cellchat, PPI.human)
}

cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

pathway <- cellchat@LR$LRsig
xlsx::write.xlsx(pathway,file=file.path('10.12/la/1',"cell_communication.xlsx"),sheetName='Sheet1',row.names=F)
write.table(pathway,file.path(project_path,paste(opt$project,'cell_communication.txt',sep='.')),sep = '\t',row.names = F,quote = F)
lapply(unique(unique(cellchat@meta$labels)),function(celltype){
  png(file.path('10.12/la/2',paste(celltype,'circle','png',sep='.')),width = 1300,height = 1300,res=300)
  print(mynetVisual_aggregate(cellchat, signaling = cellchat@netP$pathways, layout = "circle", 
                              vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7,signaling.name='ALL',
                              from=celltype))
  dev.off()
  pdf(file.path('10.12/la/2',paste(celltype,'circle','pdf',sep='.')))
  print(mynetVisual_aggregate(cellchat, signaling = cellchat@netP$pathways, layout = "circle", 
                              vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7,signaling.name='ALL',
                              from=celltype))
  dev.off()
})

pdf(file.path(project_path,paste('all_circle.pdf',sep='')))
print(netVisual_aggregate(cellchat, signaling = cellchat@netP$pathways, layout = "circle", vertex.size = groupSize,pt.title=20,vertex.label.cex = 1.7,signaling.name='ALL'))
dev.off()

#nPatterns = 5 
#cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
#cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
#png(file.path(project_path,paste('outgoing.dotplot.png',sep='')),width = 2000,height = 794,res=300)
#netAnalysis_dot(cellchat, pattern = "outgoing") + 
#  theme(plot.title = element_text(size=18),axis.text = element_text(size=15),legend.text = element_text(size=13),legend.title = element_text(size=15))
#dev.off()
##netAnalysis_dot(cellchat, pattern = "outgoing") + 
 # theme(plot.title = element_text(size=18),axis.text = element_text(size=15),legend.text = element_text(size=13),legend.title = element_text(size=15))
#dev.off()

#cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
#cellchat <- netAnalysis_signalingRole(cellchat, slot.name = "netP")
#png(file.path(project_path,paste('incoming.dotplot.png',sep='')),width = 2000,height = 794,res=300)
#netAnalysis_dot(cellchat, pattern = "incoming") + 
#  theme(plot.title = element_text(size=18),axis.text = element_text(size=15),legend.text = element_text(size=13),legend.title = element_text(size=15))
#dev.off()
#pdf(file.path(project_path,paste('incoming.dotplot.pdf',sep='')))
#netAnalysis_dot(cellchat, pattern = "incoming") + 
#  theme(plot.title = element_text(size=18),axis.text = element_text(size=15),legend.text = element_text(size=13),legend.title = element_text(size=15))
#dev.off()