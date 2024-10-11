library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(ggdendro)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)
library(concaveman)
library(scCustomize)
library(reshape)
library(dendextend)
library(biomaRt)
library(MatrixExtra)
library(rtracklayer)

##################################

##################################

expr_data <- function(data = data, genes = genes){
  
  gen_expr <- data[genes,]
  
  if (!is.vector(gen_expr)){gen_expr <- colSums(as.matrix(gen_expr))}
  
  return(gen_expr)
}

##################################


tissue_datasets <- c("Endotel_other_IL","Endotel_other_RE","Endotel_other_TC",
                     "Epithelium_IL","Epithelium_RE","Epithelium_TC",
                     "immune_IL","immune_RE","immune_TC")

cell_names_df <- as.data.frame(Reduce(rbind,sapply(tissue_datasets, function(x){
  
  filename <- list.files(path = paste0("../three_locations/",x,"/identification/"))
  filename <- filename[sapply(filename,function(x){grepl("RDS",x, fixed = T)})]
  if(length(filename)==2){filename <- filename[2]}
  readRDS(paste0("../three_locations/",x,"/identification/",filename))
  
})))

seurat_object <- readRDS("../seurat_aggregated_filtered.RDS")

############################_TISSUES_DENDROGRAMM_###############################

mean_embeddings <- t(sapply(unique(cell_names_df$cell_names), function(x){
  
  cells <- rownames(cell_names_df[cell_names_df$cell_names==x,])
  embeddings <- seurat_object@reductions$harmony@cell.embeddings[
    rownames(seurat_object@reductions$harmony@cell.embeddings)%in%cells,]
  colMeans(embeddings)
  
}))

dist_mean_embeddings <- dist(mean_embeddings, "euclidean")
clustering_mean_embeddings <- hclust(dist_mean_embeddings, method = "complete")

dendr_tissues <- dendro_data(clustering_mean_embeddings, type="rectangle") 

color_info <- sapply(label(dendr_tissues)$label, function(x){
  
  if (!is.na(as.numeric(tail(strsplit(x, "_")[[1]], n = 1)))){"other"} else 
    {sub("IL_","",sub("RE_","",sub("TC_","",x)))}
  
})

color_info <- as.factor(color_info)
levels(color_info) <- scCustomize_Palette(num_groups = length(levels(color_info)), 
                                        ggplot_default_colors = T)
levels(color_info)[levels(color_info)==names(which.max(table(color_info)))] <- "#000000"
color_info_names <- names(color_info)
color_info <- as.character(color_info)
names(color_info) <- color_info_names

dendr_tissues$labels$label <- sapply(label(dendr_tissues)$label, function(x){
  paste0(x,"::",sum(cell_names_df$cell_names==x),"::",
         length(unique(sapply(rownames(cell_names_df[cell_names_df$cell_names==x,]), 
                              function(y){
                                
                                strsplit(y,"-")[[1]][3]
                                
                                }))),"/57")
})

dendrogram_tissues <- ggplot() + 
  geom_segment(data=segment(dendr_tissues), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr_tissues), aes(x=x, y=y, label=label, 
                                   hjust=0), colour = color_info, size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        legend.position = "none")

pdf(file="dendrogram_tissues.pdf", width = 15, height = 25)

plot(dendrogram_tissues)

dev.off()


cell_names_df <- cell_names_df[match(names(seurat_object$orig.ident),rownames(cell_names_df)),]
all(names(seurat_object$orig.ident)==rownames(cell_names_df))

seurat_object$cell_names <- cell_names_df$cell_names

all_cols <- scCustomize_Palette(num_groups = length(unique(seurat_object$cell_names)), 
                               ggplot_default_colors = FALSE)

all_tissues <- DimPlot(seurat_object,
                       reduction = "umap",
                       group.by = c("cell_names"), 
                       cols = all_cols, label = T,
                       label.size = 3)

pdf(file="all_UMAP.pdf", width = 30, height = 12)

plot(all_tissues)

dev.off()

############################_CLUSTERS_DENDROGRAMM_##############################

mean_embeddings_cluster <- t(sapply(unique(cell_names_df$clusters), function(x){
  
  cells <- rownames(cell_names_df[cell_names_df$clusters==x,])
  embeddings <- seurat_object@reductions$harmony@cell.embeddings[
    rownames(seurat_object@reductions$harmony@cell.embeddings)%in%cells,]
  colMeans(embeddings)
  
}))

dist_mean_embeddings_cluster <- dist(mean_embeddings_cluster, "euclidean")
clustering_mean_embeddings_cluster <- hclust(dist_mean_embeddings_cluster, 
                                             method = "complete")

dendr_clusters <- dendro_data(clustering_mean_embeddings_cluster, type="rectangle") 


color_info_clusters <- dendr_clusters$labels$label
names(color_info_clusters) <- color_info_clusters

color_info_clusters[] <- sapply(color_info_clusters, function(x){
  
  unname(color_info[names(color_info)==unique(cell_names_df$cell_names[cell_names_df$clusters==x])])
  
})

dendr_clusters$labels$label <- sapply(label(dendr_clusters)$label, function(x){
  paste0(x,"::",sum(cell_names_df$clusters==x),"::",
         length(unique(sapply(rownames(cell_names_df[cell_names_df$clusters==x,]), 
                              function(y){
                                
                                strsplit(y,"-")[[1]][3]
                                
                              }))),"/57")
})

dendrogram_clusters <- ggplot() + 
  geom_segment(data=segment(dendr_clusters), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=label(dendr_clusters), aes(x=x, y=y, label=label, hjust=0),
            colour = color_info_clusters, size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank(),
        legend.position = "none")

pdf(file="dendrogram_clusters.pdf", width = 15, height = 35)

plot(dendrogram_clusters)

dev.off()

dendr_clusters_stat <- as.data.frame(t(sapply(label(dendr_clusters)$label,function(x){
  
  strsplit(x,"::")[[1]][2:3]
  
})))

dendr_clusters_stat[,2] <- sapply(dendr_clusters_stat[,2],function(x){
  
  strsplit(x,"/", fixed = T)[[1]][1]
  
})

dendr_clusters_stat[,] <- apply(dendr_clusters_stat,2,as.numeric)
colnames(dendr_clusters_stat) <- c("N_cells","N_samples")
dendr_clusters_stat$cells_per_sample <- dendr_clusters_stat$N_cells/dendr_clusters_stat$N_samples

summary(dendr_clusters_stat$cells_per_sample)
summary(dendr_clusters_stat$N_cells)
ggplot(dendr_clusters_stat, aes(x = cells_per_sample))+
  theme_bw()+geom_density()+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15))

nnodes(clustering_mean_embeddings_cluster)
nleaves(clustering_mean_embeddings_cluster)

################################################################################
##############################_MAKE_DATASETS_###################################
################################################################################

all(rownames(cell_names_df)==names(seurat_object$orig.ident))
seurat_object$eqtl_leaves <- cell_names_df$clusters


Genes_info <- readGFF("genes.gtf")
Genes_info <- Genes_info[Genes_info$type=="gene",
                         c("gene_id","gene_name","seqid","start","strand")]

Genes_info_filtered <- Genes_info[Genes_info$seqid%in%c(1:22, "X"),]
Genes_info_filtered$seqid <- as.character(Genes_info_filtered$seqid)
Genes_info_filtered$seqid[Genes_info_filtered$seqid=="X"] <- "23"
Genes_info_filtered$seqid <- as.numeric(Genes_info_filtered$seqid)
Genes_info_filtered <- Genes_info_filtered[order(Genes_info_filtered$seqid,
                                                 Genes_info_filtered$start),]

################################################################################

bed_list <- vector(mode = "list", length = length(clustering_mean_embeddings_cluster$labels))

for (i in 1:length(bed_list)){
  
  cells_of_cluster <- cell_names_df[cell_names_df$clusters==clustering_mean_embeddings_cluster$labels[i],]
  cells_of_donor <- as.data.frame(cbind(rownames(cells_of_cluster),
                                        sapply(rownames(cells_of_cluster), function(y){
                                          strsplit(y,"-",fixed = T)[[1]][3]
                                        })))
  colnames(cells_of_donor) <- c("cell","donor")
  cells_per_donor <- table(cells_of_donor$donor)
  
  percent_LP <- seurat_object@meta.data[rownames(cells_of_cluster),c("orig.ident","MULTI_ID")]
  percent_LP$MULTI_ID <- sapply(percent_LP$MULTI_ID, function(x){grepl("LP",x)})
  percent_LP <- as.matrix(table(percent_LP))
  
  if (dim(percent_LP)[2]==1){ if (colnames(percent_LP)=="TRUE") {
    
    cells_per_donor <- rbind(cells_per_donor,rep(1,length(cells_per_donor)))
    rownames(cells_per_donor) <- c("N_cells","percent_LP")
    
  } else {
    
    cells_per_donor <- rbind(cells_per_donor,rep(0,length(cells_per_donor)))
    rownames(cells_per_donor) <- c("N_cells","percent_LP")
    
  } } else {
    
    if (nrow(percent_LP)==1){
      
      percent_LP <- percent_LP/sum(percent_LP[,1:2])
      cells_per_donor <- rbind(cells_per_donor,percent_LP[,2])
      rownames(cells_per_donor) <- c("N_cells","percent_LP")
      
    } else {
      
      percent_LP <- percent_LP/rowSums(percent_LP[,1:2])
      cells_per_donor <- rbind(cells_per_donor,percent_LP[,2])
      rownames(cells_per_donor) <- c("N_cells","percent_LP")
      
    }
  }
  
  
  expression_of_cluster <- seurat_object@assays$RNA@counts[,cells_of_donor$cell]
  rownames(expression_of_cluster) <- unname(sapply(rownames(expression_of_cluster), function(y){
    strsplit(y,"-",fixed = T)[[1]][1]
  }))
  
  expression_of_cluster <- t_shallow(expression_of_cluster)
  expression_of_cluster <- as.data.frame(t(as.matrix(fac2sparse(cells_of_donor$donor)%*%
                                                       expression_of_cluster)))
  
  shared_genes <- intersect(rownames(expression_of_cluster),
                            Genes_info_filtered$gene_id)
  
  Genes_info_step <- Genes_info_filtered[Genes_info_filtered$gene_id%in%shared_genes,]
  expression_of_cluster <- expression_of_cluster[shared_genes,,drop = F]
  expression_of_cluster <- expression_of_cluster[match(Genes_info_step$gene_id,
                                                       rownames(expression_of_cluster)),,drop = F]
  
  bed_file <- cbind(Genes_info_step$seqid,
                    Genes_info_step$start,
                    Genes_info_step$start+1,
                    Genes_info_step$gene_id,
                    Genes_info_step$gene_id,
                    Genes_info_step$strand,
                    expression_of_cluster)
  
  colnames(bed_file)[1:6] <- c("chr","start","end","phenotype_id",
                               "phenotype_group_id","strand")
  bed_file$chr[bed_file$chr==23] <- "X"
  colnames(bed_file)[1] <- "#chr"
  
  bed_list[[i]] <- list(bed_file,cells_per_donor)
  
  
}

names(bed_list) <- clustering_mean_embeddings_cluster$labels


saveRDS(bed_list,"bed_list.rds")

dendro_mean_embeddings_cluster <- as.dendrogram(clustering_mean_embeddings_cluster)



################################################################################

# remove seurat object as we dont need it anymore
rm(seurat_object); gc()

all_levels <- partition_leaves(dendro_mean_embeddings_cluster)
saveRDS(all_levels,"all_levels.rds")

bed_list_all <- vector(mode = "list", length = length(all_levels))

for (i in 1:length(all_levels)){
  
  if (length(all_levels[[i]])==1){
    
    bed_list_all[[i]] <- bed_list[[all_levels[[i]]]]
    
  } else {
    
    ## make a bed file
    
    shared_columns <- unique(unlist(lapply(bed_list[all_levels[[i]]], function(x){
      
      colnames(x[[1]])[7:length(colnames(x[[1]]))]
      
    })))
    
    expression_mat <- sapply(shared_columns, function(x){
      
      column_step <- lapply(bed_list[all_levels[[i]]], function(y){y[[1]][,colnames(y[[1]])%in%x]})
      ifexist <- lapply(column_step, length)
      column_step <- column_step[ifexist>0]
      Reduce("+", column_step)
      
    })
    
    bed_file <- cbind(bed_list[[1]][[1]][,1:6],as.data.frame(expression_mat))
    
    ## make a cell number vector
    
    cells_per_donor <- sapply(shared_columns, function(x){
      
      column_step <- lapply(bed_list[all_levels[[i]]], function(y){y[[2]][,colnames(y[[2]])%in%x]})
      ifexist <- lapply(column_step, length)
      column_step <- column_step[ifexist>0]
      column_step <- (Reduce("cbind", column_step))
      
      if (is.matrix(column_step)){
        
        c(sum(column_step[1,]), weighted.mean(x = column_step[2,], w = column_step[1,]))
        
      } else {unname(column_step)}
      
    })
    
    rownames(cells_per_donor) <- c("N_cells","percent_LP")
    
    bed_list_all[[i]] <- list(bed_file,cells_per_donor)

  }
  
}

saveRDS(bed_list_all,"bed_list_all.rds")
saveRDS(dendro_mean_embeddings_cluster, "cluster_dendrogram.rds")



