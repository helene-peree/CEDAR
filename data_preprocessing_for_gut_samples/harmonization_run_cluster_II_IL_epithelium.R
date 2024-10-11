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
cell_names_ileum <- readRDS("../../identification/cell_names_ileum.RDS")

seurat_objects <- seurat_objects[unname(unlist(split(1:length(seurat_objects), 
                                                     sort(1:length(seurat_objects)%%5))[arguments]))]

seurat_objects <- lapply(seurat_objects, function(x){
  
  x <- subset(x, subset = cells_for_sure == "singlets_same_type")
  
  x$bioptates <- sapply(x$HTO_classification, function(y){
    
    sub(".1","", fixed = T, gsub(".2","", fixed = T, y))
    
  })
  
  x$bioptates <- sapply(x$bioptates, function(y){
    
    z <- unlist(strsplit(y, "_", fixed = T))
    
    if (length(z)==1){z} else if (length(unique(z)) == 1){z[1]} else{paste(z,collapse="_")}
    
  })
  
  x$locations <- sapply(x$bioptates, function(y){
    
    strsplit(y, "-", fixed = T)[[1]][1]
    
    
  })
  
  if ("IL"%in%x$locations){
    
    # filter epithelium cells
    
    x <- x[,colnames(x)%in%names(cell_names_ileum)]
    cell_names_ileum_step <- cell_names_ileum[names(cell_names_ileum)%in%names(x$orig.ident)]
    cell_names_ileum_step <- cell_names_ileum_step[match(names(x$orig.ident),
                                                         names(cell_names_ileum_step))]
    x$cell_names <- cell_names_ileum_step
    x <- subset(x, subset = cell_names == "epithelium")
    
    #
    
    x <- subset(x = x, subset = locations == "IL")
    x <- NormalizeData(x, assay = "RNA", normalization.method = "LogNormalize")
    x <- ScaleData(x)
    x <- FindVariableFeatures(x, selection.method = "mean.var.plot", nfeatures = 3000)
    
  } else {x <- 0}
  
  
})

seurat_objects <- seurat_objects[unlist(lapply(seurat_objects, function(x){typeof(x)!="double"}))]

saveRDS(seurat_objects, paste0("seurat_part_IL_epithelium_",arguments,".RDS"))






