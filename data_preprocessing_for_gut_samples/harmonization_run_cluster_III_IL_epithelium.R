library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)


seurat_objects <- readRDS("seurat_part_IL_epithelium_1.RDS")
seurat_objects <- c(seurat_objects,readRDS("seurat_part_IL_epithelium_2.RDS"))
seurat_objects <- c(seurat_objects,readRDS("seurat_part_IL_epithelium_3.RDS"))
seurat_objects <- c(seurat_objects,readRDS("seurat_part_IL_epithelium_4.RDS"))
seurat_objects <- c(seurat_objects,readRDS("seurat_part_IL_epithelium_5.RDS"))

more_49_cells <- unlist(lapply(seurat_objects, function(x){length(x$cell_names)>=50}))

features_integr <- SelectIntegrationFeatures(object.list = seurat_objects[more_49_cells], 
                                             assay = rep("RNA",
                                                         sum(more_49_cells)))
seurat_aggregated <- merge(x = seurat_objects[[1]], 
                           y = seurat_objects[2:length(seurat_objects)], merge.data=TRUE)

rm(seurat_objects); gc()

print("finish preparation")

VariableFeatures(seurat_aggregated) <- features_integr
seurat_aggregated <- ScaleData(seurat_aggregated)
seurat_aggregated <- RunPCA(object = seurat_aggregated, assay = "RNA", npcs = 50)
seurat_aggregated <- RunHarmony(object = seurat_aggregated,
                                assay.use = "RNA",
                                reduction = "pca",
                                dims.use = 1:50,
                                plot_convergence = TRUE,
                                group.by.vars = "orig.ident")

seurat_aggregated <- RunUMAP(object = seurat_aggregated, 
                             assay = "RNA", reduction = "harmony", dims = 1:50)
seurat_aggregated <- RunTSNE(seurat_aggregated, reduction = "harmony", dims = 1:50)

seurat_aggregated <- FindNeighbors(object = seurat_aggregated, assay = "RNA", 
                                   reduction = "harmony", dims = 1:50)
seurat_aggregated <- FindClusters(object = seurat_aggregated, resolution = 1.5)

saveRDS(seurat_aggregated, "seurat_aggregated_IL_epithelium.RDS")

plot_singlets_for_sure <- DimPlot(seurat_aggregated, reduction = "tsne",
                                  group.by = c("seurat_clusters","hash.ID","MULTI_ID"))

plot_singlets_for_sure_UMAP <- DimPlot(seurat_aggregated, reduction = "umap",
                                       group.by = c("seurat_clusters","hash.ID","MULTI_ID"))

###

pdf("tSNE_filtered_IL_epithelium.pdf",
    height = 5, width = 16)

plot(plot_singlets_for_sure)

dev.off()

###

pdf("UMAP_filtered_IL_epithelium.pdf",
    height = 5, width = 16)

plot(plot_singlets_for_sure_UMAP)

dev.off()

###
