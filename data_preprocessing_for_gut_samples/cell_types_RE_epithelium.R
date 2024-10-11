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


##################################

##############################_RECTUM_###########################################

seurat_object_RE <- readRDS("seurat_aggregated_RE_epithelium.RDS")

seurat_object_RE@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_RE@assays$RNA@data),
                                                               function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][2]})

seurat_object_RE@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_RE@assays$RNA@data),
                                                                    function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][1]})

seurat_object_RE@assays$RNA@meta.features$full_name <- rownames(seurat_object_RE@assays$RNA@data)
rownames(seurat_object_RE@assays$RNA@data) <- seurat_object_RE@assays$RNA@meta.features$gene_names
rownames(seurat_object_RE@assays$RNA@counts) <- seurat_object_RE@assays$RNA@meta.features$gene_names

RE_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_RE$seurat_clusters)), 
                               ggplot_default_colors = FALSE)

DimPlot(seurat_object_RE, 
        reduction = "umap",
        group.by = c("seurat_clusters"), 
        cols = RE_cols)




####____________________________LEVEL_3_____________________________________####

#############################_PANETH_CELLS_##########################################

genes_test_paneth_cells <- c("DEFA5","DEFA6","REG3A")
seurat_object_RE$paneth_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                          genes_test_paneth_cells)



#############################_GOBLET_CELLS_##########################################

genes_test_goblet_cells <- c("MUC2")
seurat_object_RE$goblet_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                          genes_test_goblet_cells)



#############################_STEM_CELLS_#########################################

genes_test_stem_cells <- c("LGR5","ASCL2","OLFM4","SMOC2",
                         "RGMB")
seurat_object_RE$stem_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                        genes_test_stem_cells)



#############################_Enteroendocrine_CELLS_############################

genes_test_enteroend_cells <- c("CHGA","CHGB","NEUROD1") 
seurat_object_RE$enteroend_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                              genes_test_enteroend_cells)


###########################_BEST4_CELLS_########################################

genes_test_BEST4_cells <- c("BEST4","CA7","OTOP2") 
seurat_object_RE$BEST4_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                         genes_test_BEST4_cells)



################################_TUFT_CELLS_#####################################

genes_test_tuft_cells <- c("POU2F3","IRAG2","TRPM5","SH2D6") 
seurat_object_RE$tuft_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                        genes_test_tuft_cells)


################################_TA_CELLS_#####################################

genes_test_TA_cells <- c("MKI67","TOP2A","PCNA") 
seurat_object_RE$TA_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                      genes_test_TA_cells)



###########################_IL_entero_CELLS_#########################################

genes_test_IL_entero_cells <- c("FABP2","FABP6") 
seurat_object_RE$IL_entero_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                             genes_test_IL_entero_cells)

###########################_RE_entero_CELLS_#########################################

genes_test_RE_entero_cells <- c("B3GNT7","LCN2","PRAC1") 
seurat_object_RE$RE_entero_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                             genes_test_RE_entero_cells)

###########################_TC_entero_CELLS_#########################################

genes_test_TC_entero_cells <- c("TXN","CD24","CA2") 
seurat_object_RE$TC_entero_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                             genes_test_TC_entero_cells)

###########################_TCRE_entero_CELLS_#########################################

genes_test_TCRE_entero_cells <- c("TMEM37","CD177","CA4") #CA2 excluded
seurat_object_RE$TCRE_entero_vect <- expr_data(seurat_object_RE@assays$RNA@data, 
                                             genes_test_TCRE_entero_cells)

##_dotplot

RE_epithelium_3_lvl_dotplot <- DotPlot(seurat_object_RE, 
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

pdf(file="identification/RE_epithelium_3_lvl_dotplot.pdf", width = 24, height = 12)

plot(RE_epithelium_3_lvl_dotplot)

dev.off()

##_cells

RE_epithelium_3_lvl_plot <- DimPlot(seurat_object_RE, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = RE_cols, label = T) +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "paneth_vect") +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "goblet_vect") +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "stem_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "enteroend_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "BEST4_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "tuft_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "TA_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "IL_entero_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "RE_entero_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "TC_entero_vect") + 
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "TCRE_entero_vect") + 
  plot_layout(ncol = 3)

##

pdf(file="identification/RE_epithelium_3_lvl_plot.pdf", width = 16, height = 20)

plot(RE_epithelium_3_lvl_plot)

dev.off()

##

cell_names_RE_epithelial <- as.numeric(as.character(seurat_object_RE$seurat_clusters))
names(cell_names_RE_epithelial) <- names(seurat_object_RE$seurat_clusters)

cell_names_RE_epithelial <- sapply(cell_names_RE_epithelial, function(x){
  
  if (x%in%c(30,10)) {x <- "RE_BEST4"} else 
    if (x==29) {x <- "RE_Tuft"} else 
      if (x%in%c(15)) {x <- "RE_TA"} else
        if (x%in%c(0,4,18,21)) {x <- "RE_enterocytes_TCRE"} else 
          if (x%in%c(3,5,26,32)) {x <- "RE_goblet"} else
            if (x%in%c(28)) {x <- "RE_enteroendocrine"} else
              if (x%in%c(15,6,22,1,8,11,34,19,14,24,9,10,26,13,27,2,17,23,20,16)) {x <- "RE_enterocytes_RE"} else
                if (x%in%c(12)) {x <- "RE_enterocytes_IL"} else
                  if (x%in%c(7)) {x <- "RE_stem"} else
                         {x <- paste0("RE_epithelial_",x)}
  
})


seurat_object_RE$cell_subnames <- cell_names_RE_epithelial

pdf(file="identification/RE_epithelial_3_lvl_tissues.pdf", width = 9, height = 8)

DimPlot(seurat_object_RE, 
        reduction = "umap",
        group.by = c("cell_subnames"), 
        cols = RE_cols, label = T,
        label.size = 2) 

dev.off()

cell_names_RE_epithelial_table <- cbind(as.numeric(as.character(seurat_object_RE$seurat_clusters)),
                                        cell_names_RE_epithelial)
colnames(cell_names_RE_epithelial_table) <- c("clusters","cell_names")

cell_names_RE_epithelial_table[,1] <- paste0("RE_epithelial_",
                                             cell_names_RE_epithelial_table[,1])

saveRDS(cell_names_RE_epithelial_table,"identification/cell_names_RE_epithelial.RDS")

################################################################################

seurat_object_RE <- RunUMAP(object = seurat_object_RE, n.components = 3L,
                            assay = "RNA", reduction = "harmony", dims = 1:50)

##########################################
