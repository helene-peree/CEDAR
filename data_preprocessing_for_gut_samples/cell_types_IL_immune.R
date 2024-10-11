library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)
library(concaveman)
library(scCustomize)
library(reshape)

library(rgl)
library(magick)

##################################

##################################

expr_data <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){gen_expr <- colSums(as.matrix(gen_expr))}
  
  return(gen_expr)
}

expr_data_NK_T <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){
    
    gen_expr_more_0 <- apply(gen_expr,c(1,2),function(x){x>0})
    gen_expr_more_0 <- apply(gen_expr_more_0, 2, all)
    gen_expr <- colSums(as.matrix(gen_expr))
    gen_expr <- gen_expr*gen_expr_more_0
    
    }
  
  return(gen_expr)
}


##################################

##############################_ILEUM_###########################################

seurat_object_IL <- readRDS("seurat_aggregated_IL_immune.RDS")

seurat_object_IL@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_IL@assays$RNA@data),
                                                               function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][2]})

seurat_object_IL@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_IL@assays$RNA@data),
                                                                    function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][1]})

seurat_object_IL@assays$RNA@meta.features$full_name <- rownames(seurat_object_IL@assays$RNA@data)
rownames(seurat_object_IL@assays$RNA@data) <- seurat_object_IL@assays$RNA@meta.features$gene_names
rownames(seurat_object_IL@assays$RNA@counts) <- seurat_object_IL@assays$RNA@meta.features$gene_names

IL_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_IL$seurat_clusters)), 
                               ggplot_default_colors = FALSE)

DimPlot(seurat_object_IL, 
        reduction = "umap",
        group.by = c("seurat_clusters"), 
        cols = IL_cols)


####____________________________LEVEL_3_____________________________________####

#############################_T_CELLS_##########################################

genes_test_Tcells <- c("CD3G","CD3D","CD3E")
seurat_object_IL$Tcells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_Tcells)

genes_test_Tcells.CD4 <- c("CD4")
seurat_object_IL$Tcells_vect.CD4  <- expr_data(seurat_object_IL@assays$RNA@data, 
                                               genes_test_Tcells.CD4)

genes_test_Tcells.CD8 <- c("CD8A","CD8B")
seurat_object_IL$Tcells_vect.CD8  <- expr_data(seurat_object_IL@assays$RNA@data, 
                                               genes_test_Tcells.CD8)

genes_test_T_cells.NK_T <- c("CD3D","GZMA","NKG7","PRF1") 
seurat_object_IL$T_cells_vect.NK_T <- expr_data_NK_T(seurat_object_IL@assays$RNA@data, 
                                                genes_test_T_cells.NK_T)

genes_test_T_cells.gdT <- c("TRGC1") 
seurat_object_IL$T_cells_vect.gdT <- expr_data(seurat_object_IL@assays$RNA@data, 
                                                     genes_test_T_cells.gdT)

genes_test_T_cells.MAIT <- c("TRAV1-1", "TRAV1-2","TRAV2","SLC4A10") 
seurat_object_IL$T_cells_vect.MAIT <- expr_data(seurat_object_IL@assays$RNA@data, 
                                                    genes_test_T_cells.MAIT)


#############################_B_CELLS_##########################################

genes_test_Bcells <- c("MS4A1","CD19","CD79A","CD79B")
seurat_object_IL$Bcells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_Bcells)

genes_test_Bcells.Naive <- c("IGHD")
seurat_object_IL$Bcells_vect.Naive <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_Bcells.Naive)

genes_test_Bcells.Memory <- c("CD27")
seurat_object_IL$Bcells_vect.Memory <- expr_data(seurat_object_IL@assays$RNA@data, 
                                                genes_test_Bcells.Memory)

genes_test_Bcells.Plasma <- c("SDC1","IGHA1","IGHA2","CD38")
seurat_object_IL$Bcells_vect.Plasma <- expr_data(seurat_object_IL@assays$RNA@data, 
                                                 genes_test_Bcells.Plasma)

#############################_NK_CELLS_#########################################

genes_test_NK_cells <- c("NKG7","EOMES","NCAM1","PRF1",
                         "KIR3DL1") ### add another KIRs?
seurat_object_IL$NK_cells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_NK_cells)


#############################_ILCs_CELLS_#######################################

genes_test_ILCs_cells_negative <- c("CD3G","CD3D","CD3E","MS4A1","CD19","CD79A","CD79B",
                           "CD14","FUT4","FCGR3A","FCGR3B") 

seurat_object_IL$ILCs_cells_vect_neg <- expr_data(seurat_object_IL@assays$RNA@counts, 
                                                  genes_test_ILCs_cells_negative)

genes_test_ILCs_cells <- c("TBX21","IKZF3","IL7R","GATA3","MAF","PTGDR2","HPGDS",
                           "KIT","IL1R1","IL23R","RORC") 

seurat_object_IL$ILCs_cells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                              genes_test_ILCs_cells)

seurat_object_IL$ILCs_cells_vect <- unlist(mapply(function(x,y){
  
  if(y>0){x <- 0}else{x <- x}
  x
  
}, x = seurat_object_IL$ILCs_cells_vect, y = seurat_object_IL$ILCs_cells_vect_neg))



###########################_MYELOIDS_CELLS_#####################################

genes_test_myel_cells <- c("LYZ","ITGAM") 
seurat_object_IL$myel_cells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                              genes_test_myel_cells)

genes_test_monocytes_cells <- c("FCN1","S100A4","S100A6") 
seurat_object_IL$myel_cells_vect.monocytes <- expr_data_NK_T(seurat_object_IL@assays$RNA@data, 
                                                        genes_test_monocytes_cells)


################################_pDC_CELLS_#####################################

genes_test_pDC_cells <- c("JCHAIN","CLEC4C") 
seurat_object_IL$pDC_cells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_pDC_cells)

################################_oDC_CELLS_#####################################

genes_test_oDC_cells <- c("CLEC9A","CLEC10A","LAMP3") 
seurat_object_IL$oDC_cells_vect <- expr_data_NK_T(seurat_object_IL@assays$RNA@data, 
                                             genes_test_oDC_cells)



###########################_MAST_CELLS_#########################################

genes_test_mast_cells <- c("GATA2","CPA3","HPGDS") 
seurat_object_IL$mast_cells_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_mast_cells)

##_dotplot

IL_immune_3_lvl_dotplot <- DotPlot(seurat_object_IL, 
                                   features = unique(c(genes_test_Tcells,
                                                       genes_test_Tcells.CD4,
                                                       genes_test_Tcells.CD8,
                                                       genes_test_T_cells.NK_T,
                                                       genes_test_T_cells.gdT,
                                                       genes_test_T_cells.MAIT,
                                                       genes_test_Bcells,
                                                       genes_test_Bcells.Naive,
                                                       genes_test_Bcells.Memory,
                                                       genes_test_Bcells.Plasma,
                                                       genes_test_NK_cells,
                                                       genes_test_myel_cells,
                                                       genes_test_monocytes_cells,
                                                       genes_test_pDC_cells,
                                                       genes_test_oDC_cells,
                                                       genes_test_mast_cells)))

pdf(file="identification/IL_immune_3_lvl_dotplot.pdf", width = 18, height = 12)

plot(IL_immune_3_lvl_dotplot)

dev.off()

##_cells

IL_immune_3_lvl_plot <- DimPlot(seurat_object_IL, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = IL_cols, label = T) +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Tcells_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Tcells_vect.CD4") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Tcells_vect.CD8") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "T_cells_vect.NK_T") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "T_cells_vect.gdT") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "T_cells_vect.MAIT") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Bcells_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Bcells_vect.Naive") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Bcells_vect.Memory") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "Bcells_vect.Plasma") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "NK_cells_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "ILCs_cells_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "myel_cells_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "myel_cells_vect.monocytes") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "pDC_cells_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "oDC_cells_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "mast_cells_vect") + 
  plot_layout(ncol = 3)

##

pdf(file="identification/IL_immune_3_lvl_plot.pdf", width = 16, height = 20)

plot(IL_immune_3_lvl_plot)

dev.off()

##

cell_names_IL_immune <- as.numeric(as.character(seurat_object_IL$seurat_clusters))
names(cell_names_IL_immune) <- names(seurat_object_IL$seurat_clusters)

cell_names_IL_immune <- sapply(cell_names_IL_immune, function(x){
  
  if (x%in%c(5,10,8,1,2,0,17,13,7,20,11,24,21)) {x <- "IL_T_cells"} else 
    if (x%in%c(4,25,22,3,27,23,6,18,15)) {x <- "IL_B_cells"} else 
      if (x%in%c()) {x <- "IL_ILC"} else
        if (x%in%c(16,9)) {x <- "IL_myeloids"} else 
          if (x%in%c(19)) {x <- "IL_pDC"} else
            if (x%in%c(12)) {x <- "IL_mast_cells"} else
              {x <- paste0("IL_immune_",x)}
  
})


seurat_object_IL$cell_subnames <- cell_names_IL_immune

pdf(file="identification/IL_immune_3_lvl_tissues.pdf", width = 9, height = 8)

DimPlot(seurat_object_IL, 
        reduction = "umap",
        group.by = c("cell_subnames"), 
        cols = IL_cols, label = T) 

dev.off()

cell_names_IL_immune_table <- cbind(as.numeric(as.character(seurat_object_IL$seurat_clusters)),
                                    cell_names_IL_immune)
colnames(cell_names_IL_immune_table) <- c("clusters","cell_names")

cell_names_IL_immune_table[,1] <- paste0("IL_immune_",
                                         cell_names_IL_immune_table[,1])

saveRDS(cell_names_IL_immune_table,"identification/cell_names_IL_immune.RDS")

################_DEEPER

cell_names_IL_immune_4 <- as.numeric(as.character(seurat_object_IL$seurat_clusters))
names(cell_names_IL_immune_4) <- names(seurat_object_IL$seurat_clusters)

cell_names_IL_immune_4 <- sapply(cell_names_IL_immune_4, function(x){
  
  if (x%in%c(8,7,3)) {x <- "IL_T_cells.CD4"} else 
    if (x%in%c(0,13,6)) {x <- "IL_T_cells.CD8"} else 
      if (x%in%c(5)) {x <- "IL_T_cells.MAIT"} else
        if (x%in%c(1)) {x <- "IL_B_cells.Naive"} else 
          if (x%in%c(20,2,17,4,16,14)) {x <- "IL_B_cells.Memory"} else
            if (x%in%c(15)) {x <- "IL_myeloids.monocytes"} else
            {x <- paste0("IL_immune_",x)}
  
})


seurat_object_IL$cell_subnames_4 <- cell_names_IL_immune_4

pdf(file="identification/IL_immune_4_lvl_tissues.pdf", width = 9, height = 8)

DimPlot(seurat_object_IL, 
        reduction = "umap",
        group.by = c("cell_subnames_4"), 
        cols = IL_cols, label = T) 

dev.off()

cell_names_IL_immune_4_table <- cbind(as.numeric(as.character(seurat_object_IL$seurat_clusters)),
                                      cell_names_IL_immune_4)
colnames(cell_names_IL_immune_4_table) <- c("clusters","cell_names")

cell_names_IL_immune_4_table[,1] <- paste0("IL_immune_",
                                         cell_names_IL_immune_4_table[,1])

saveRDS(cell_names_IL_immune_4_table,"identification/cell_names_IL_immune_4.RDS")


###############################################################################


seurat_object_IL <- RunUMAP(object = seurat_object_IL, n.components = 3L,
                            assay = "RNA", reduction = "harmony", dims = 1:50)

##########################################











