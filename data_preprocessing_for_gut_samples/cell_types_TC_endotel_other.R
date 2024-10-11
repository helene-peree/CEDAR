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

##################################

##################################

expr_data <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){gen_expr <- colSums(as.matrix(gen_expr))}
  
  return(gen_expr)
}


##################################

##############################_TCEUM_###########################################

seurat_object_TC <- readRDS("seurat_aggregated_TC_endotel_other.RDS")

seurat_object_TC@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                               function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][2]})

seurat_object_TC@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                                    function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][1]})

seurat_object_TC@assays$RNA@meta.features$full_name <- rownames(seurat_object_TC@assays$RNA@data)
rownames(seurat_object_TC@assays$RNA@data) <- seurat_object_TC@assays$RNA@meta.features$gene_names

TC_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_TC$seurat_clusters)), 
                               ggplot_default_colors = FALSE)

DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("seurat_clusters"), 
        cols = TC_cols)




####____________________________LEVEL_3_____________________________________####

#############################_FIBROBLASTS_CELLS_##########################################

genes_test_fibroblasts <- c("DCN")
seurat_object_TC$fibroblasts_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                               genes_test_fibroblasts)

#############################_smooth_m_CELLS_##########################################

genes_test_smooth_m <- c("DES","CNN1")
seurat_object_TC$smooth_m_vect <- expr_data(seurat_object_TC@assays$RNA@data,
                                            genes_test_smooth_m)


#############################_Myofibroblast_CELLS_##########################################

genes_test_myofibroblast <- c("ACTA2","TAGLN")
seurat_object_TC$myofibroblast_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                                 genes_test_myofibroblast)

#############################_Enteric_glia_CELLS_#########################################

genes_test_enteric_glia <- c("FOXD3","MPZ","PLP1") 
seurat_object_TC$enteric_glia_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                                genes_test_enteric_glia)

#############################_Pericytes_CELLS_#######################################

genes_test_pericytes_cells <- c("RGS5","ACTG2","NOTCH3","MCAM") 
seurat_object_TC$pericytes_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                             genes_test_pericytes_cells)

#############################_Endoteliocytes_CELLS_#############################

genes_test_endoteliocytes_cells <- c("VWF","PECAM1","CDH5") 
seurat_object_TC$endoteliocytes_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                                  genes_test_endoteliocytes_cells)

#############################_Neural_CELLS_#####################################

genes_test_neural_cells <- c("ETV1","BNC2") 
seurat_object_TC$neural_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                          genes_test_neural_cells)


#############################_mesodermal_CELLS_######################################

genes_test_mesodermal_cells <- c("HAND1", "HAND2", "PITX2", "ZEB2")
seurat_object_TC$mesodermal_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                              genes_test_mesodermal_cells)

#############################_stromal_CELLS_######################################

genes_test_stromal_cells <- c("ADAMDEC1", "PDGFRA", "BMP4", "C7",
                              "MMP1", "MMP3", "PDPN", "COL7A1", "CHI3L1")
seurat_object_TC$stromal_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                           genes_test_stromal_cells)

#############################_Cajal_CELLS_######################################

genes_test_cajal_cells <- c("KIT", "ANO1")
seurat_object_TC$cajal_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                         genes_test_cajal_cells)

#############################_Mesothelial_CELLS_#################################

genes_test_mesotelial_cells <- c("KRT19", "LRRN4", "UPK3B")
seurat_object_TC$mesotelial_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                              genes_test_mesotelial_cells)

#############################_Immune_contam_CELLS_##############################

genes_test_immune_contam_cells <- c("PTPRC")

seurat_object_TC$immune_contam_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                                 genes_test_immune_contam_cells)

#############################_Epithelial_contam_CELLS_##############################

genes_test_epi_contam_cells <- c("EPCAM")

seurat_object_TC$epi_contam_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                              genes_test_epi_contam_cells)



##_dotplot

TC_endotel_other_3_lvl_dotplot <- DotPlot(seurat_object_TC, 
                                          features = unique(c(genes_test_fibroblasts,
                                                              genes_test_smooth_m,
                                                              genes_test_myofibroblast,
                                                              genes_test_enteric_glia,
                                                              genes_test_pericytes_cells,
                                                              genes_test_endoteliocytes_cells,
                                                              genes_test_neural_cells,
                                                              genes_test_mesodermal_cells,
                                                              genes_test_stromal_cells,
                                                              genes_test_cajal_cells,
                                                              genes_test_mesotelial_cells,
                                                              genes_test_immune_contam_cells,
                                                              genes_test_epi_contam_cells)))

pdf(file="identification/TC_endotel_other_3_lvl_dotplot.pdf", width = 24, height = 12)

plot(TC_endotel_other_3_lvl_dotplot)

dev.off()

##_cells

TC_endotel_other_3_lvl_plot <- DimPlot(seurat_object_TC, 
                                       reduction = "umap",
                                       group.by = c("seurat_clusters"), 
                                       cols = TC_cols, label = T) +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "fibroblasts_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "smooth_m_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "myofibroblast_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "enteric_glia_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "pericytes_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "endoteliocytes_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "neural_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "mesodermal_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "stromal_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "cajal_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "mesotelial_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "immune_contam_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "epi_contam_vect") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/TC_endotel_other_3_lvl_plot.pdf", width = 10, height = 24)

plot(TC_endotel_other_3_lvl_plot)

dev.off()

##

cell_names_TC_endotel_other <- as.numeric(as.character(seurat_object_TC$seurat_clusters))
names(cell_names_TC_endotel_other) <- names(seurat_object_TC$seurat_clusters)

cell_names_TC_endotel_other <- sapply(cell_names_TC_endotel_other, function(x){
  
  if (x==17) {x <- "TC_enteric_glia"} else 
    if (x%in%c(14,22)) {x <- "TC_endoteliocytes"} else
    {x <- paste0("TC_endotel_other_",x)}
  
})

seurat_object_TC$cell_subnames <- cell_names_TC_endotel_other

pdf(file="identification/TC_endotel_other_3_lvl_tissues.pdf", width = 9, height = 8)

DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("cell_subnames"), 
        cols = TC_cols, label = T,
        label.size = 2) 

dev.off()

cell_names_TC_endotel_other_table <- cbind(as.numeric(as.character(seurat_object_TC$seurat_clusters)),
                                           cell_names_TC_endotel_other)
colnames(cell_names_TC_endotel_other_table) <- c("clusters","cell_names")

cell_names_TC_endotel_other_table[,1] <- paste0("TC_endotel_other_",
                                                cell_names_TC_endotel_other_table[,1])

saveRDS(cell_names_TC_endotel_other_table,"identification/cell_names_TC_endotel_other.RDS")



