arguments <- as.numeric(commandArgs(trailingOnly = T))

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)

seurat_objects <- readRDS("seurat_objects_all.RDS")

seurat_objects <- seurat_objects[unname(unlist(split(1:length(seurat_objects), 
                                                     sort(1:length(seurat_objects)%%5))[arguments]))]

seurat_objects <- lapply(seurat_objects, function(x){
  
  x <- subset(x, subset = cells_for_sure == "singlets_same_type")
  x <- NormalizeData(x, assay = "RNA", normalization.method = "LogNormalize")
  x <- ScaleData(x)
  x <- FindVariableFeatures(x, selection.method = "mean.var.plot", nfeatures = 3000)
  
  
})

saveRDS(seurat_objects, paste0("seurat_part_",arguments,".RDS"))
