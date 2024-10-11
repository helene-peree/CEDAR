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

expr_data <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){gen_expr <- colSums(as.matrix(gen_expr))}
  
  return(gen_expr)
}


cell_names_ileum <- readRDS("cell_names_ileum.RDS")
cell_names_ileum <- cell_names_ileum[cell_names_ileum=="epithelium"]
cell_names_ileum[] <- paste0(cell_names_ileum,"_IL")

cell_names_rectum <- readRDS("cell_names_rectum.RDS")
cell_names_rectum <- cell_names_rectum[cell_names_rectum=="epithelium"]
cell_names_rectum[] <- paste0(cell_names_rectum,"_RE")

cell_names_transversum <- readRDS("cell_names_transversum.RDS")
cell_names_transversum <- cell_names_transversum[cell_names_transversum=="epithelium"]
cell_names_transversum[] <- paste0(cell_names_transversum,"_TC")

cell_names_epithelium <- c(cell_names_ileum,cell_names_rectum,cell_names_transversum)

saveRDS(cell_names_epithelium,"cell_names_epithelium.RDS")

