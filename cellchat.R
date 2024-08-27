library(CellChat)
library(patchwork)
library(Seurat)
data.input = sce.all@assays$RNA@data
meta.data =  sce.all@meta.data
table(meta.data$seurat_annotations)
head(meta.data)
cellchat <- createCellChat(object = data.input, 
                           meta = meta.data, 
                           group.by = "seurat_annotations")
cellchat <- addMeta(cellchat, meta = meta.data)
levels(cellchat@idents) # show factor levels of the cell labels
cellchat <- setIdent(cellchat, ident.use = "seurat_annotations") 
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group
levels(cellchat@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
future::plan('multiprocess', workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat) 
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)