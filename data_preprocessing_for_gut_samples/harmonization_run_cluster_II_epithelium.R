arguments <- as.numeric(commandArgs(trailingOnly = T))

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)

seurat_objects <- readRDS("../../seurat_objects_all.RDS")
cell_names_epithelium <- readRDS("cell_names_epithelium.RDS")

seurat_objects <- seurat_objects[unname(unlist(split(1:length(seurat_objects), 
                                                     sort(1:length(seurat_objects)%%5))[arguments]))]

seurat_objects <- lapply(seurat_objects, function(x){
  
  x <- x[,colnames(x)%in%names(cell_names_epithelium)]
  x <- NormalizeData(x, assay = "RNA", normalization.method = "LogNormalize")
  x <- ScaleData(x)
  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", nfeatures = 3000)

})

seurat_objects <- seurat_objects[unlist(lapply(seurat_objects, function(x){typeof(x)!="double"}))]

saveRDS(seurat_objects, paste0("seurat_part_epithelium_",arguments,".RDS"))






