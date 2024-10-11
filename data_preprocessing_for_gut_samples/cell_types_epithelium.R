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

expr_data <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){gen_expr <- colSums(as.matrix(gen_expr))}
  
  return(gen_expr)
}


##################################

##############################_ILEUM_###########################################

seurat_object <- readRDS("seurat_aggregated_epithelium.RDS")

cell_names_epithelium <- readRDS("cell_names_epithelium.RDS")
cell_names_epithelium <- cell_names_epithelium[match(names(seurat_object$orig.ident),
                                                     names(cell_names_epithelium))]

all(names(seurat_object$orig.ident)==names(cell_names_epithelium))
seurat_object$cell_names <- cell_names_epithelium

IL_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object$seurat_clusters)), 
                               ggplot_default_colors = FALSE)

DimPlot(seurat_object, 
        reduction = "umap",
        group.by = c("cell_names"))




####____________________________LEVEL_3_____________________________________####

#############################_PANETH_CELLS_##########################################

genes_test_paneth_cells <- c("DEFA5","DEFA6","REG3A")
seurat_object$paneth_vect <- expr_data(seurat_object@assays$RNA@data, 
                                          genes_test_paneth_cells)



#############################_GOBLET_CELLS_##########################################

genes_test_goblet_cells <- c("MUC2")
seurat_object$goblet_vect <- expr_data(seurat_object@assays$RNA@data, 
                                          genes_test_goblet_cells)



#############################_STEM_CELLS_#########################################

genes_test_stem_cells <- c("LGR5","ASCL2","OLFM4","SMOC2",
                           "RGMB")
seurat_object$stem_vect <- expr_data(seurat_object@assays$RNA@data, 
                                        genes_test_stem_cells)



#############################_Enteroendocrine_CELLS_############################

genes_test_enteroend_cells <- c("CHGA","CHGB","NEUROD1") 
seurat_object$enteroend_vect <- expr_data(seurat_object@assays$RNA@data, 
                                             genes_test_enteroend_cells)


###########################_BEST4_CELLS_########################################

genes_test_BEST4_cells <- c("BEST4","CA7","OTOP2") 
seurat_object$BEST4_vect <- expr_data(seurat_object@assays$RNA@data, 
                                         genes_test_BEST4_cells)



################################_TUFT_CELLS_#####################################

genes_test_tuft_cells <- c("POU2F3","IRAG2","TRPM5","SH2D6") 
seurat_object$tuft_vect <- expr_data(seurat_object@assays$RNA@data, 
                                        genes_test_tuft_cells)


################################_TA_CELLS_#####################################

genes_test_TA_cells <- c("MKI67","TOP2A","PCNA") 
seurat_object$TA_vect <- expr_data(seurat_object@assays$RNA@data, 
                                      genes_test_TA_cells)



###########################_IL_entero_CELLS_#########################################

genes_test_IL_entero_cells <- c("FABP2","FABP6") 
seurat_object$IL_entero_vect <- expr_data(seurat_object@assays$RNA@data, 
                                             genes_test_IL_entero_cells)

###########################_RE_entero_CELLS_#########################################

genes_test_RE_entero_cells <- c("B3GNT7","LCN2","PRAC1") 
seurat_object$RE_entero_vect <- expr_data(seurat_object@assays$RNA@data, 
                                             genes_test_RE_entero_cells)

###########################_TC_entero_CELLS_#########################################

genes_test_TC_entero_cells <- c("TXN","CD24","CA2") 
seurat_object$TC_entero_vect <- expr_data(seurat_object@assays$RNA@data, 
                                             genes_test_TC_entero_cells)

###########################_TCRE_entero_CELLS_#########################################

genes_test_TCRE_entero_cells <- c("TMEM37","CD177","CA4") #CA2 excluded
seurat_object$TCRE_entero_vect <- expr_data(seurat_object@assays$RNA@data, 
                                               genes_test_TCRE_entero_cells)

###########################_DSC_CELLS_#########################################

genes_test_DSC_cells <- c("REG4") #CA2 excluded
seurat_object$DSC_vect <- expr_data(seurat_object@assays$RNA@data, 
                                    genes_test_DSC_cells)

##_dotplot

IL_epithelium_3_lvl_dotplot <- DotPlot(seurat_object, 
                                       features = unique(c(genes_test_paneth_cells,
                                                           genes_test_goblet_cells,
                                                           genes_test_stem_cells,
                                                           genes_test_enteroend_cells,
                                                           genes_test_BEST4_cells,
                                                           genes_test_tuft_cells,
                                                           genes_test_TA_cells,
                                                           genes_test_IL_entero_cells,
                                                           genes_test_RE_entero_cells,
                                                           genes_test_TC_entero_cells,
                                                           genes_test_TCRE_entero_cells,
                                                           genes_test_DSC_cells)))

pdf(file="identification/epithelium_3_lvl_dotplot.pdf", width = 24, height = 12)

plot(IL_epithelium_3_lvl_dotplot)

dev.off()

##_cells

IL_epithelium_3_lvl_plot <- DimPlot(seurat_object, 
                                    reduction = "umap",
                                    group.by = c("seurat_clusters"), 
                                    cols = IL_cols) +
  DimPlot(seurat_object, reduction = "umap", group.by = c("cell_names")) +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "paneth_vect") +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "goblet_vect") +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "stem_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "enteroend_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "BEST4_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "tuft_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "TA_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "IL_entero_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "RE_entero_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "TC_entero_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "TCRE_entero_vect") + 
  FeaturePlot(seurat_object, reduction = "umap",
              features = "DSC_vect") + 
  plot_layout(ncol = 3)

##

pdf(file="identification/epithelium_3_lvl_plot.pdf", width = 16, height = 20)

plot(IL_epithelium_3_lvl_plot)

dev.off()

##



seurat_object <- RunUMAP(object = seurat_object, n.components = 3L,
                            assay = "RNA", reduction = "harmony", dims = 1:50)


