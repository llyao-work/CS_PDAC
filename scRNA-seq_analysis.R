library(tidyverse)
library(Seurat)

myfile <- list.files('Seurat1/')
myfile[1:6]
mysample <- unlist(lapply(myfile,function(x){strsplit(x,'_')[[1]][2]}))
mysample[1:6]
fileformat <-  unlist(lapply(myfile,function(x){strsplit(x,'_')[[1]][3]}))
fileformat[1:6]
for(i in 1: length(myfile)){
  dir.create(paste0('Seurat1/',mysample[i]))
  file.copy(paste0('Seurat1/',myfile[i]),
            paste0('Seurat1/',mysample[i],'/',fileformat[i]))
  unlink(paste0('Seurat1/',myfile[i]))
}
samples1=list.files('Seurat1/')
samples1
sceList1 = lapply(samples1,function(pro){
  folder=file.path("Seurat1/",pro)
  CreateSeuratObject(counts = Read10X(folder),
                     project = pro )
})
View(sceList1)
sce.all=merge(x=sceList1[[1]],
              y=sceList1[ -1 ],
              add.cell.ids = samples1)
as.data.frame(sce.all@assays$RNA@counts[1:10, 1:2])

sce <- ace.all
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce)
sce <- RunPCA(sce,assay="RNA")
library(harmony)
sce <- RunHarmony(sce,group.by.vars="orig.ident" )
sce <- RunTSNE(sce, reduction= "harmony")
sce <- RunUMAP(sce, reduction = "harmony")
sce <- FindNeighbors(sce, reduction = "harmony") 
sce <- FindClusters(object = sce, verbose = T, resolution = 0.5)
sce@meta.data$seurat_clusters <- sce@meta.data$RNA_snn_res.0.5