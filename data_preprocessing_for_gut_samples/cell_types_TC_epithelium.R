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

##############################_TRANSVERSUM_###########################################

seurat_object_TC <- readRDS("seurat_aggregated_TC_epithelium.RDS")

seurat_object_TC@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                               function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][2]})

seurat_object_TC@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                                    function(x){strsplit(sub("-","*",x),"*",fixed = T)[[1]][1]})

seurat_object_TC@assays$RNA@meta.features$full_name <- rownames(seurat_object_TC@assays$RNA@data)
rownames(seurat_object_TC@assays$RNA@data) <- seurat_object_TC@assays$RNA@meta.features$gene_names
rownames(seurat_object_TC@assays$RNA@counts) <- seurat_object_TC@assays$RNA@meta.features$gene_names

TC_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_TC$seurat_clusters)), 
                               ggplot_default_colors = FALSE)

DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("seurat_clusters"), 
        cols = TC_cols)




####____________________________LEVEL_3_____________________________________####

#############################_PANETH_CELLS_##########################################

genes_test_paneth_cells <- c("DEFA5","DEFA6","REG3A")
seurat_object_TC$paneth_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                          genes_test_paneth_cells)



#############################_GOBLET_CELLS_##########################################

genes_test_goblet_cells <- c("MUC2")
seurat_object_TC$goblet_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                          genes_test_goblet_cells)



#############################_STEM_CELLS_#########################################

genes_test_stem_cells <- c("LGR5","ASCL2","OLFM4","SMOC2",
                         "RGMB")
seurat_object_TC$stem_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                        genes_test_stem_cells)



#############################_Enteroendocrine_CELLS_############################

genes_test_enteroend_cells <- c("CHGA","CHGB","NEUROD1") 
seurat_object_TC$enteroend_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                              genes_test_enteroend_cells)


###########################_BEST4_CELLS_########################################

genes_test_BEST4_cells <- c("BEST4","CA7","OTOP2") 
seurat_object_TC$BEST4_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                         genes_test_BEST4_cells)



################################_TUFT_CELLS_#####################################

genes_test_tuft_cells <- c("POU2F3","IRAG2","TRPM5","SH2D6") 
seurat_object_TC$tuft_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                        genes_test_tuft_cells)


################################_TA_CELLS_#####################################

genes_test_TA_cells <- c("MKI67","TOP2A","PCNA") 
seurat_object_TC$TA_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                      genes_test_TA_cells)



###########################_IL_entero_CELLS_#########################################

genes_test_IL_entero_cells <- c("FABP2","FABP6") 
seurat_object_TC$IL_entero_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                             genes_test_IL_entero_cells)

###########################_RE_entero_CELLS_#########################################

genes_test_RE_entero_cells <- c("B3GNT7","LCN2","PRAC1") 
seurat_object_TC$RE_entero_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                             genes_test_RE_entero_cells)

###########################_TC_entero_CELLS_#########################################

genes_test_TC_entero_cells <- c("TXN","CD24","CA2") 
seurat_object_TC$TC_entero_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                             genes_test_TC_entero_cells)

###########################_TCRE_entero_CELLS_#########################################

genes_test_TCRE_entero_cells <- c("TMEM37","CD177","CA4") #CA2 excluded
seurat_object_TC$TCRE_entero_vect <- expr_data(seurat_object_TC@assays$RNA@data, 
                                             genes_test_TCRE_entero_cells)

##_dotplot

TC_epithelium_3_lvl_dotplot <- DotPlot(seurat_object_TC, 
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

pdf(file="identification/TC_epithelium_3_lvl_dotplot.pdf", width = 24, height = 12)

plot(TC_epithelium_3_lvl_dotplot)

dev.off()

##_cells

TC_epithelium_3_lvl_plot <- DimPlot(seurat_object_TC, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = TC_cols, label = T) +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "paneth_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "goblet_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "stem_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "enteroend_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "BEST4_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "tuft_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "TA_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "IL_entero_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "RE_entero_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "TC_entero_vect") + 
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "TCRE_entero_vect") + 
  plot_layout(ncol = 3)

##

pdf(file="identification/TC_epithelium_3_lvl_plot.pdf", width = 16, height = 20)

plot(TC_epithelium_3_lvl_plot)

dev.off()

##

cell_names_TC_epithelial <- as.numeric(as.character(seurat_object_TC$seurat_clusters))
names(cell_names_TC_epithelial) <- names(seurat_object_TC$seurat_clusters)

cell_names_TC_epithelial <- sapply(cell_names_TC_epithelial, function(x){
  
  if (x%in%c(17,21,36)) {x <- "TC_BEST4"} else 
    if (x==22) {x <- "TC_Tuft"} else 
      if (x%in%c(13,5)) {x <- "TC_TA"} else
        if (x%in%c(0,9,1,23)) {x <- "TC_enterocytes_TCRE"} else 
          if (x%in%c(20,2,8,4)) {x <- "TC_goblet"} else
            if (x%in%c(30)) {x <- "TC_enteroendocrine"} else
              if (x%in%c(12,10,32,29,11,15,3,7,28,35,18,6,16,27,34)) {x <- "TC_enterocytes_TC"} else
                if (x%in%c(19,25,24)) {x <- "TC_enterocytes_IL"} else
                  if (x%in%c(14)) {x <- "TC_stem"} else
                    {x <- paste0("TC_epithelial_",x)}
  
})


seurat_object_TC$cell_subnames <- cell_names_TC_epithelial

pdf(file="identification/TC_epithelial_3_lvl_tissues.pdf", width = 9, height = 8)

plot(DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("cell_subnames"), 
        cols = TC_cols, label = T,
        label.size = 2))

dev.off()

cell_names_TC_epithelial_table <- cbind(as.numeric(as.character(seurat_object_TC$seurat_clusters)),
                                        cell_names_TC_epithelial)
colnames(cell_names_TC_epithelial_table) <- c("clusters","cell_names")

cell_names_TC_epithelial_table[,1] <- paste0("TC_epithelial_",
                                             cell_names_TC_epithelial_table[,1])

saveRDS(cell_names_TC_epithelial_table,"identification/cell_names_TC_epithelial.RDS")

#######################################################################


seurat_object_TC <- RunUMAP(object = seurat_object_TC, n.components = 3L,
                            assay = "RNA", reduction = "harmony", dims = 1:50)

##########################################



