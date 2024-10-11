library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)
library(concaveman)
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


##################################

##############################_ILEUM_###########################################

seurat_object_IL <- readRDS("seurat_aggregated_IL_epithelium.RDS")

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

#############################_PANETH_CELLS_##########################################

genes_test_paneth_cells <- c("DEFA5","DEFA6","REG3A")
seurat_object_IL$paneth_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_paneth_cells)



#############################_GOBLET_CELLS_##########################################

genes_test_goblet_cells <- c("MUC2")
seurat_object_IL$goblet_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                          genes_test_goblet_cells)



#############################_STEM_CELLS_#########################################

genes_test_stem_cells <- c("LGR5","ASCL2","OLFM4","SMOC2",
                         "RGMB")
seurat_object_IL$stem_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                        genes_test_stem_cells)



#############################_Enteroendocrine_CELLS_############################

genes_test_enteroend_cells <- c("CHGA","CHGB","NEUROD1") 
seurat_object_IL$enteroend_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                              genes_test_enteroend_cells)


###########################_BEST4_CELLS_########################################

genes_test_BEST4_cells <- c("BEST4","CA7","OTOP2") 
seurat_object_IL$BEST4_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                         genes_test_BEST4_cells)



################################_TUFT_CELLS_#####################################

genes_test_tuft_cells <- c("POU2F3","IRAG2","TRPM5","SH2D6") 
seurat_object_IL$tuft_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                        genes_test_tuft_cells)


################################_TA_CELLS_#####################################

genes_test_TA_cells <- c("MKI67","TOP2A","PCNA") 
seurat_object_IL$TA_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                      genes_test_TA_cells)



###########################_IL_entero_CELLS_#########################################

genes_test_IL_entero_cells <- c("FABP2","FABP6") 
seurat_object_IL$IL_entero_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_IL_entero_cells)

###########################_RE_entero_CELLS_#########################################

genes_test_RE_entero_cells <- c("B3GNT7","LCN2","PRAC1") 
seurat_object_IL$RE_entero_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_RE_entero_cells)

###########################_TC_entero_CELLS_#########################################

genes_test_TC_entero_cells <- c("TXN","CD24","CA2") 
seurat_object_IL$TC_entero_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_TC_entero_cells)

###########################_TCRE_entero_CELLS_#########################################

genes_test_TCRE_entero_cells <- c("TMEM37","CD177","CA4") #CA2 excluded
seurat_object_IL$TCRE_entero_vect <- expr_data(seurat_object_IL@assays$RNA@data, 
                                             genes_test_TCRE_entero_cells)

##_dotplot

IL_epithelium_3_lvl_dotplot <- DotPlot(seurat_object_IL, 
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
                                                           genes_test_TCRE_entero_cells)))

pdf(file="identification/IL_epithelium_3_lvl_dotplot.pdf", width = 24, height = 12)

plot(IL_epithelium_3_lvl_dotplot)

dev.off()

##_cells

IL_epithelium_3_lvl_plot <- DimPlot(seurat_object_IL, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = IL_cols, label = T) +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "paneth_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "goblet_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "stem_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "enteroend_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "BEST4_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "tuft_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "TA_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "IL_entero_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "RE_entero_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "TC_entero_vect") + 
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "TCRE_entero_vect") + 
  plot_layout(ncol = 4)

##

pdf(file="identification/IL_epithelium_3_lvl_plot.pdf", width = 20, height = 16)

plot(IL_epithelium_3_lvl_plot)

dev.off()

##

##

jpeg(file="identification/IL_epithelium_3_lvl_plot.jpeg", width = 1800, height = 1200)

plot(IL_epithelium_3_lvl_plot)

dev.off()

##

cell_names_IL_epithelial <- as.numeric(as.character(seurat_object_IL$seurat_clusters))
names(cell_names_IL_epithelial) <- names(seurat_object_IL$seurat_clusters)

cell_names_IL_epithelial <- sapply(cell_names_IL_epithelial, function(x){
  
  if (x==23) {x <- "IL_BEST4"} else 
    if (x==15) {x <- "IL_Tuft"} else 
      if (x==8) {x <- "IL_TA"} else
        if (x==19) {x <- "IL_enterocytes_TCRE"} else 
            if (x%in%c(6,20,34,18,13)) {x <- "IL_goblet"} else
              if (x%in%c(28,32)) {x <- "IL_enteroendocrine"} else
                if (x%in%c(0,2,22,11,17,26,1,5,3,
                           36,12,31,30,21,4,7,25)) {x <- "IL_enterocytes_IL"} else
                             if (x==27){x <- "IL_paneth"} else 
                               if (x==24){x <- "IL_enterocytes_RE"} else
                                 if (x==9){x <- "IL_stem"} else {x <- paste0("IL_epithelial_",x)}
   
})


seurat_object_IL$cell_subnames <- cell_names_IL_epithelial

pdf(file="identification/IL_epithelial_3_lvl_tissues.pdf", width = 9, height = 8)

DimPlot(seurat_object_IL, 
        reduction = "umap",
        group.by = c("cell_subnames"), 
        cols = IL_cols, label = T,
        label.size = 2) 

dev.off()

cell_names_IL_epithelial_table <- cbind(as.numeric(as.character(seurat_object_IL$seurat_clusters)),
                                           cell_names_IL_epithelial)
colnames(cell_names_IL_epithelial_table) <- c("clusters","cell_names")

cell_names_IL_epithelial_table[,1] <- paste0("IL_epithelial_",
                                             cell_names_IL_epithelial_table[,1])

saveRDS(cell_names_IL_epithelial_table,"identification/cell_names_IL_epithelial.RDS")

################################################################################

seurat_object_IL <- RunUMAP(object = seurat_object_IL, n.components = 3L,
                            assay = "RNA", reduction = "harmony", dims = 1:50)

seurat_object_IL

all.markers <- FindAllMarkers(object = seurat_object_IL,
                              assay = "RNA")

##########################################

