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


#########################_MAIN_3_CELLTYPES_#####################################

##############################_ILEUM_###########################################

seurat_object_IL <- readRDS("three_locations/seurat_aggregated_IL_filtered.RDS")

seurat_object_IL@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_IL@assays$RNA@data),
                                                function(x){strsplit(x,"-",fixed = T)[[1]][2]})

seurat_object_IL@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_IL@assays$RNA@data),
                                                               function(x){strsplit(x,"-",fixed = T)[[1]][1]})

seurat_object_IL@assays$RNA@meta.features$full_name <- rownames(seurat_object_IL@assays$RNA@data)
rownames(seurat_object_IL@assays$RNA@data) <- seurat_object_IL@assays$RNA@meta.features$gene_names

genes <- c("PTPRC","EPCAM","VWF")

genes_test <- "PTPRC"
seurat_object_IL$PTPRC_vect <- expr_data(seurat_object_IL@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object_IL$EPCAM_vect <- expr_data(seurat_object_IL@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_IL$VWF_vect <- expr_data(seurat_object_IL@assays$RNA@data, genes_test)

##

##

IL_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_IL$seurat_clusters)), 
                               ggplot_default_colors = FALSE)
IL_3_types_plot <- DimPlot(seurat_object_IL, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = IL_cols, label = T) + ggtitle('Seurat_clusters') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "PTPRC_vect") + ggtitle('PTPRC') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "EPCAM_vect") + ggtitle('EPCAM') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "VWF_vect") + ggtitle("VWF, PECAM1, CDH5") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/IL_3_types_plot.pdf", width = 12, height = 12)

plot(IL_3_types_plot)

dev.off()

##

IL_3_types_Dotplot <- DotPlot(seurat_object_IL, assay = "RNA", features = genes)+RotatedAxis()

##

pdf(file="identification/IL_3_types_dotplot.pdf", width = 8, height = 12)

plot(IL_3_types_Dotplot)

dev.off()

#######################_TISSUE_IDENTIFICATION_IL_###############################

cell_names_ileum <- as.numeric(as.character(seurat_object_IL$seurat_clusters))
names(cell_names_ileum) <- names(seurat_object_IL$seurat_clusters)

PTPRC_clusters <- c(1,6,11,26,46,29,15,36,23,44,42,4,39,19,25,13,50)
EPCAM_clusters <- c(0,2,51,3,22,7,37,49,12,8,9,47,17,35,18,52,21,31,16,10,28,31,27,32,53,38,24,20,43,12)
VWF_other_cluster <- c(30,5,45,33,14,40,41,48,34) #23 and 8 moved to immune

cell_names_ileum <- sapply(cell_names_ileum, function(x){
  
  if (x%in%PTPRC_clusters) {x <- "immune"} else 
    if (x%in%EPCAM_clusters) {x <- "epithelium"} else
      if (x%in%VWF_other_cluster){x <- "endotel_other"} else {x}
  
})

seurat_object_IL$cell_names_ileum <- cell_names_ileum

# seurat_object_IL <- ScaleData(seurat_object_IL,
#                               features = rownames(seurat_object_IL@assays$RNA@data))

genes_test <- "PTPRC"
seurat_object_IL$PTPRC_vect <- expr_data(seurat_object_IL@assays$RNA@scale.data, genes_test)

genes_test <- "EPCAM"
seurat_object_IL$EPCAM_vect <- expr_data(seurat_object_IL@assays$RNA@scale.data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_IL$VWF_vect <- expr_data(seurat_object_IL@assays$RNA@scale.data, genes_test)


##

IL_3_types_plot <- DimPlot(seurat_object_IL, 
                           reduction = "umap",
                           group.by = c("cell_names_ileum"), label = T) +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "PTPRC_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "EPCAM_vect") +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "VWF_vect") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/IL_3_tissues_plot.pdf", width = 12, height = 12)

plot(IL_3_types_plot)

dev.off()

##

saveRDS(cell_names_ileum,"identification/cell_names_ileum.RDS")


################################################################################

IL_3_types_cluster_expr <- reshape(IL_3_types_Dotplot$data[,2:4], idvar = "id",
                                   timevar = "features.plot", direction = "wide")

IL_3_types_cluster_expr$simple_clusters <- apply(IL_3_types_cluster_expr[,2:4], 
                                                 1, function(x){
  
  max_cluster <- which.max(x)
  
  if (x[max_cluster]<=20){"strange"} else 
    if (any(x[-which.max(x)]>=20)){"mixed"}else{names(max_cluster)}
  
})

##

transcript_counts_IL <- as.matrix(seurat_object_IL@assays$RNA@counts[genes,])
transcript_counts_IL_bool <- apply(transcript_counts_IL,2,function(x){x>0})
transcript_counts_IL_double_pos <- transcript_counts_IL[,colSums(transcript_counts_IL_bool)>1]
transcript_counts_IL_negative <- transcript_counts_IL[,colSums(transcript_counts_IL_bool)==0]

dim(transcript_counts_IL)[2]
dim(transcript_counts_IL_double_pos)[2]
dim(transcript_counts_IL_negative)[2]

min_IL_double_pos <- unlist(apply(transcript_counts_IL_double_pos,2,
                           function(x){min(x[x>0])}))

max_IL_double_pos <- unlist(apply(transcript_counts_IL_double_pos,2,
                                  function(x){max(x[x>0])}))

summary(min_IL_double_pos)
plot(density(min_IL_double_pos))

summary(max_IL_double_pos)
plot(density(max_IL_double_pos))

seurat_object_IL$doublePositive <- as.factor(colSums(transcript_counts_IL_bool)>1)
seurat_object_IL$nCount_RNA_log <- log(seurat_object_IL$nCount_RNA, 10)

DimPlot(seurat_object_IL, 
        reduction = "umap",
        group.by = c("doublePositive")) + ggtitle('doublePositive') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "nCount_RNA_log") + ggtitle('nCount_RNA_log') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "nFeature_RNA") + ggtitle('nFeature_RNA') +
  FeaturePlot(seurat_object_IL, reduction = "umap",
              features = "percent.mt") + ggtitle('percent.mt') + 
  plot_layout(ncol = 2)





##############################_RECTUM_###########################################

seurat_object_RE <- readRDS("three_locations/seurat_aggregated_RE_filtered.RDS")

seurat_object_RE@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_RE@assays$RNA@data),
                                                               function(x){strsplit(x,"-",fixed = T)[[1]][2]})

seurat_object_RE@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_RE@assays$RNA@data),
                                                                    function(x){strsplit(x,"-",fixed = T)[[1]][1]})

seurat_object_RE@assays$RNA@meta.features$full_name <- rownames(seurat_object_RE@assays$RNA@data)
rownames(seurat_object_RE@assays$RNA@data) <- seurat_object_RE@assays$RNA@meta.features$gene_names

genes <- c("PTPRC","EPCAM","VWF")


genes_test <- "PTPRC"
seurat_object_RE$PTPRC_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object_RE$EPCAM_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_RE$VWF_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)

##

RE_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_RE$seurat_clusters)), 
                               ggplot_default_colors = FALSE)
RE_3_types_plot <- DimPlot(seurat_object_RE, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = RE_cols, label = T) + ggtitle('Seurat_clusters') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "PTPRC_vect") + ggtitle('PTPRC') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "EPCAM_vect") + ggtitle('EPCAM') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "VWF_vect") + ggtitle("VWF, PECAM1, CDH5") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/RE_3_types_plot.pdf", width = 12, height = 12)

plot(RE_3_types_plot)

dev.off()

##

RE_3_types_Dotplot <- DotPlot(seurat_object_RE, features = genes)+RotatedAxis()

##

pdf(file="identification/RE_3_types_dotplot.pdf", width = 8, height = 12)

plot(RE_3_types_Dotplot)

dev.off()

#######################_TISSUE_IDENTIFICATION_RE_###############################

cell_names_rectum <- as.numeric(as.character(seurat_object_RE$seurat_clusters))
names(cell_names_rectum) <- names(seurat_object_RE$seurat_clusters)

PTPRC_clusters <- c(38,35,8,32,13,19,18,36,40,28,39)
EPCAM_clusters <- c(0,12,21,4,7,1,2,17,22,11,16,6,5,9,14,26,29,11,31,23,34,42,37,15,3,27)
VWF_other_cluster <- c(30,24,25,33,20,10,41) 

cell_names_rectum <- sapply(cell_names_rectum, function(x){
  
  if (x%in%PTPRC_clusters) {x <- "immune"} else 
    if (x%in%EPCAM_clusters) {x <- "epithelium"} else
      if (x%in%VWF_other_cluster){x <- "endotel_other"} else {x}
  
})

seurat_object_RE$cell_names_rectum <- cell_names_rectum

# seurat_object_RE <- ScaleData(seurat_object_RE,
#                               features = rownames(seurat_object_RE@assays$RNA@data))

genes_test <- "PTPRC"
seurat_object_RE$PTPRC_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object_RE$EPCAM_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_RE$VWF_vect <- expr_data(seurat_object_RE@assays$RNA@data, genes_test)


##

RE_3_types_plot <- DimPlot(seurat_object_RE, 
                           reduction = "umap",
                           group.by = c("cell_names_rectum"), label = T) +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "PTPRC_vect") +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "EPCAM_vect") +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "VWF_vect") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/RE_3_tissues_plot.pdf", width = 12, height = 12)

plot(RE_3_types_plot)

dev.off()

##

saveRDS(cell_names_rectum,"identification/cell_names_rectum.RDS")


################################################################################

RE_3_types_cluster_expr <- reshape(RE_3_types_Dotplot$data[,2:4], idvar = "id",
                                   timevar = "features.plot", direction = "wide")

RE_3_types_cluster_expr$simple_clusters <- apply(RE_3_types_cluster_expr[,2:4], 
                                                 1, function(x){
                                                   
                                                   max_cluster <- which.max(x)
                                                   
                                                   if (x[max_cluster]<=20){"strange"} else 
                                                     if (any(x[-which.max(x)]>=20)){"mixed"}else{names(max_cluster)}
                                                   
                                                 })

##

transcript_counts_RE <- as.matrix(seurat_object_RE@assays$RNA@counts[genes,])
transcript_counts_RE_bool <- apply(transcript_counts_RE,2,function(x){x>0})
transcript_counts_RE_double_pos <- transcript_counts_RE[,colSums(transcript_counts_RE_bool)>1]
transcript_counts_RE_negative <- transcript_counts_RE[,colSums(transcript_counts_RE_bool)==0]

dim(transcript_counts_RE)[2]
dim(transcript_counts_RE_double_pos)[2]
dim(transcript_counts_RE_negative)[2]

min_RE_double_pos <- unlist(apply(transcript_counts_RE_double_pos,2,
                                  function(x){min(x[x>0])}))

max_RE_double_pos <- unlist(apply(transcript_counts_RE_double_pos,2,
                                  function(x){max(x[x>0])}))

summary(min_RE_double_pos)
plot(density(min_RE_double_pos))

summary(max_RE_double_pos)
plot(density(max_RE_double_pos))

seurat_object_RE$doublePositive <- as.factor(colSums(transcript_counts_RE_bool)>1)
seurat_object_RE$nCount_RNA_log <- log(seurat_object_RE$nCount_RNA, 10)

DimPlot(seurat_object_RE, 
        reduction = "umap",
        group.by = c("doublePositive")) + ggtitle('doublePositive') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "nCount_RNA_log") + ggtitle('nCount_RNA_log') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "nFeature_RNA") + ggtitle('nFeature_RNA') +
  FeaturePlot(seurat_object_RE, reduction = "umap",
              features = "percent.mt") + ggtitle('percent.mt') + 
  plot_layout(ncol = 2)

################################################################################



##############################_TRANSVERSUM_#####################################

seurat_object_TC <- readRDS("three_locations/seurat_aggregated_TC_filtered.RDS")

seurat_object_TC@assays$RNA@meta.features$gene_names <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                               function(x){strsplit(x,"-",fixed = T)[[1]][2]})

seurat_object_TC@assays$RNA@meta.features$ensembl_gene_ID <- sapply(rownames(seurat_object_TC@assays$RNA@data),
                                                                    function(x){strsplit(x,"-",fixed = T)[[1]][1]})

seurat_object_TC@assays$RNA@meta.features$full_name <- rownames(seurat_object_TC@assays$RNA@data)
rownames(seurat_object_TC@assays$RNA@data) <- seurat_object_TC@assays$RNA@meta.features$gene_names

genes <- c("PTPRC","EPCAM","VWF")

genes_test <- "PTPRC"
seurat_object_TC$PTPRC_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object_TC$EPCAM_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_TC$VWF_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)

##

TC_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object_TC$seurat_clusters)), 
                               ggplot_default_colors = FALSE)
TC_3_types_plot <- DimPlot(seurat_object_TC, 
                           reduction = "umap",
                           group.by = c("seurat_clusters"), 
                           cols = TC_cols, label = T) + ggtitle('Seurat_clusters') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "PTPRC_vect") + ggtitle('PTPRC') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "EPCAM_vect") + ggtitle('EPCAM') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "VWF_vect") + ggtitle("VWF, PECAM1, CDH5") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/TC_3_types_plot.pdf", width = 12, height = 12)

plot(TC_3_types_plot)

dev.off()

##

TC_3_types_Dotplot <- DotPlot(seurat_object_TC, features = genes)+RotatedAxis()

##

pdf(file="identification/TC_3_types_dotplot.pdf", width = 8, height = 12)

plot(TC_3_types_Dotplot)

dev.off()

#######################_TISSUE_IDENTIFICATION_TC_###############################

cell_names_transversum <- as.numeric(as.character(seurat_object_TC$seurat_clusters))
names(cell_names_transversum) <- names(seurat_object_TC$seurat_clusters)

PTPRC_clusters <- c(42,15,18,5,17,32,44,12,28,23,41)
EPCAM_clusters <- c(27,22,41,9,2,21,4,14,19,3,0,25,16,36,24,13,1,48,34,33,26,43,45,8,31)
VWF_other_cluster <- c(29,44,12,28,38,20,30,11,7,10,6,40,23,35,37,39) 

cell_names_transversum <- sapply(cell_names_transversum, function(x){
  
  if (x%in%PTPRC_clusters) {x <- "immune"} else 
    if (x%in%EPCAM_clusters) {x <- "epithelium"} else
      if (x%in%VWF_other_cluster){x <- "endotel_other"} else {x}
  
})

seurat_object_TC$cell_names_transversum <- cell_names_transversum

# seurat_object_TC <- ScaleData(seurat_object_TC,
#                               features = rownames(seurat_object_TC@assays$RNA@data))

genes_test <- "PTPRC"
seurat_object_TC$PTPRC_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object_TC$EPCAM_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object_TC$VWF_vect <- expr_data(seurat_object_TC@assays$RNA@data, genes_test)


##

TC_3_types_plot <- DimPlot(seurat_object_TC, 
                           reduction = "umap",
                           group.by = c("cell_names_transversum"), label = T) +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "PTPRC_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "EPCAM_vect") +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "VWF_vect") + 
  plot_layout(ncol = 2)

##

pdf(file="identification/TC_3_tissues_plot.pdf", width = 12, height = 12)

plot(TC_3_types_plot)

dev.off()

##

saveRDS(cell_names_transversum,"identification/cell_names_transversum.RDS")


################################################################################


TC_3_types_cluster_expr <- reshape(TC_3_types_Dotplot$data[,2:4], idvar = "id",
                                   timevar = "features.plot", direction = "wide")

TC_3_types_cluster_expr$simple_clusters <- apply(TC_3_types_cluster_expr[,2:4], 
                                                 1, function(x){
                                                   
                                                   max_cluster <- which.max(x)
                                                   
                                                   if (x[max_cluster]<=20){"strange"} else 
                                                     if (any(x[-which.max(x)]>=20)){"mixed"}else{names(max_cluster)}
                                                   
                                                 })

##

transcript_counts_TC <- as.matrix(seurat_object_TC@assays$RNA@counts[genes,])
transcript_counts_TC_bool <- apply(transcript_counts_TC,2,function(x){x>0})
transcript_counts_TC_double_pos <- transcript_counts_TC[,colSums(transcript_counts_TC_bool)>1]
transcript_counts_TC_negative <- transcript_counts_TC[,colSums(transcript_counts_TC_bool)==0]

dim(transcript_counts_TC)[2]
dim(transcript_counts_TC_double_pos)[2]
dim(transcript_counts_TC_negative)[2]

min_TC_double_pos <- unlist(apply(transcript_counts_TC_double_pos,2,
                                  function(x){min(x[x>0])}))

max_TC_double_pos <- unlist(apply(transcript_counts_TC_double_pos,2,
                                  function(x){max(x[x>0])}))

summary(min_TC_double_pos)
plot(density(min_TC_double_pos))

summary(max_TC_double_pos)
plot(density(max_TC_double_pos))

################################################################################

seurat_object_TC$doublePositive <- as.factor(colSums(transcript_counts_TC_bool)>1)
seurat_object_TC$nCount_RNA_log <- log(seurat_object_TC$nCount_RNA, 10)

DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("doublePositive"), 
        cols = TC_cols) + ggtitle('doublePositive') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "nCount_RNA_log") + ggtitle('nCount_RNA_log') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "nFeature_RNA") + ggtitle('nFeature_RNA') +
  FeaturePlot(seurat_object_TC, reduction = "umap",
              features = "percent.mt") + ggtitle('percent.mt') + 
  plot_layout(ncol = 2)


seurat_object_TC <- FindClusters(object = seurat_object_TC, resolution = 0.5)

DimPlot(seurat_object_TC, 
        reduction = "umap",
        group.by = c("seurat_clusters"), 
        cols = TC_cols) + ggtitle('seurat_clusters') 

markers_TC <- FindAllMarkers(seurat_object_TC)

################################################################################

cell_names_ileum <- readRDS("identification/cell_names_ileum.RDS")
cell_names_rectum <- readRDS("identification/cell_names_rectum.RDS")
cell_names_transversum <- readRDS("identification/cell_names_transversum.RDS")

cell_names_ileum[] <- paste0(cell_names_ileum,"_IL")
cell_names_rectum[] <- paste0(cell_names_rectum,"_RE")
cell_names_transversum[] <- paste0(cell_names_transversum,"_TC")

##################################


seurat_object <- readRDS("seurat_aggregated_filtered.RDS")

seurat_object$bioptates <- sapply(seurat_object$HTO_classification, function(x){

  sub(".1","", fixed = T, gsub(".2","", fixed = T, x))

})

seurat_object$bioptates <- sapply(seurat_object$bioptates, function(x){

  y <- unlist(strsplit(x, "_", fixed = T))

  if (length(y)==1){y} else if (length(unique(y)) == 1){y[1]} else{paste(y,collapse="_")}

})

seurat_object$locations <- sapply(seurat_object$bioptates, function(x){

  strsplit(x, "-", fixed = T)[[1]][1]


})

genes_test <- "PTPRC"
seurat_object$PTPRC_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object$EPCAM_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object$VWF_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

##

cell_names <- c(cell_names_ileum,cell_names_rectum,cell_names_transversum)
cell_names <- cell_names[match(names(seurat_object$orig.ident),names(cell_names))]

seurat_object$cell_names <- cell_names

##

all_3_types_plot <- DimPlot(seurat_object,
                            reduction = "umap",
                            group.by = c("cell_names"), label = T) +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "PTPRC_vect") +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "EPCAM_vect") +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "VWF_vect") + 
  plot_layout(ncol = 2)

pdf(file="identification/ALL_3_tissues_plot.pdf", width = 12, height = 12)

plot(all_3_types_plot)

dev.off()



bioptat_types <- names(summary(as.factor(seurat_object$bioptates)))

#######################_overclassification_#####################################

seurat_object <- FindClusters(object = seurat_object, resolution = 1.5)

################################################################################

genes <- c("PTPRC","EPCAM","VWF")

genes_test <- "PTPRC"
seurat_object$PTPRC_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

genes_test <- "EPCAM"
seurat_object$EPCAM_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

genes_test <- c("VWF","PECAM1","CDH5")
seurat_object$VWF_vect <- expr_data(seurat_object@assays$RNA@data, genes_test)

##

all_cols <- scCustomize_Palette(num_groups = length(levels(seurat_object$seurat_clusters)),
                               ggplot_default_colors = FALSE)

all_3_types_plot <- DimPlot(seurat_object,
                           reduction = "umap",
                           group.by = c("locations"),
                           cols = all_cols, label = T) + ggtitle('locations') +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "PTPRC_vect") + ggtitle('PTPRC') +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "EPCAM_vect") + ggtitle('EPCAM') +
  FeaturePlot(seurat_object, reduction = "umap",
              features = "VWF_vect") + ggtitle("VWF, PECAM1, CDH5") +
  plot_layout(ncol = 2)

##

pdf(file="identification/all_3_types_plot.pdf", width = 12, height = 12)

plot(all_3_types_plot)

dev.off()

##

all_3_types_Dotplot <- DotPlot(seurat_object, features = genes)+RotatedAxis()

##

pdf(file="identification/all_3_types_dotplot.pdf", width = 8, height = 12)

plot(all_3_types_Dotplot)

dev.off()

##

all_3_types_cluster_expr <- reshape(all_3_types_Dotplot$data[,2:4], idvar = "id",
                                   timevar = "features.plot", direction = "wide")

all_3_types_cluster_expr$simple_clusters <- apply(all_3_types_cluster_expr[,2:4],
                                                 1, function(x){

                                                   max_cluster <- which.max(x)

                                                   if (x[max_cluster]<=20){"strange"} else
                                                     if (any(x[-which.max(x)]>=20)){"mixed"}else{names(max_cluster)}

                                                 })

##

cell_statistics <- readRDS("cell_statistics.RDS")

write.csv(cell_statistics,"cell_statistics.csv", quote = F)

cell_statistics_summary <- as.data.frame(t(t(apply(cell_statistics[,-1], 2, mean))))
colnames(cell_statistics_summary) <- "mean"
cell_statistics_summary$sum <- apply(cell_statistics[,-1], 2, sum)

write.csv(cell_statistics_summary,
          "cell_statistics_summary.csv", quote = F)





