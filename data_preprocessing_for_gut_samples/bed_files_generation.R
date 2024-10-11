library(dendextend)
library(compositions)
library(ggplot2)
library(factoextra)
library(ape)
library(ggdendro)



all_datasets <- readRDS("../cell_types_dendrogramm/bed_list_all.rds")
cluster_dendro <- readRDS("../cell_types_dendrogramm/cluster_dendrogram.rds")
structure_dendro <- readRDS("../cell_types_dendrogramm/all_levels.rds")

structure_dendro[unlist(lapply(structure_dendro, function(x){any(x%in%c("IL_epithelial_27"))}))]


dataset_type <- unlist(lapply(structure_dendro, function(x){
  
  if (length(x)>1){paste0("Node_len_", length(x), "_")} else {paste0("leaf_",x,"_")}
  
}))

names(structure_dendro) <- paste0(dataset_type,"pos_",c(1:length(structure_dendro)))
names(all_datasets) <- names(structure_dendro)

##################_MAKE_DENDROGRAM_WITH_DATASETS_NAMES_#########################

cluster_dendro_d <- dendro_data(cluster_dendro, type = "rectangle")
cluster_dendro_nodes <- as.data.frame(get_nodes_xy(cluster_dendro))

## the object with the information 
## about which leafs are included to a node
structure_dendro
##

cluster_dendro_nodes$Names <- names(structure_dendro) 


cluster_dendro_d$labels$label <- sapply(cluster_dendro_d$labels$label,function(x){""})

pdf(file="dendrogram_Names.pdf", width = 15, height = 35)

plot(ggplot() + 
       geom_segment(data=segment(cluster_dendro_d), aes(x=x, y=y, xend=xend, yend=yend)) + 
       geom_text(data=label(cluster_dendro_d), aes(x=x, y=y, label=label, hjust=0),
                 size=3) +
       geom_text(data=cluster_dendro_nodes, aes(x=V1, y=V2, label=Names, hjust=0),
                 size=3, color = "red") +
       coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
       theme(axis.line.y=element_blank(),
             axis.ticks.y=element_blank(),
             axis.text.y=element_blank(),
             axis.title.y=element_blank(),
             panel.background=element_rect(fill="white"),
             panel.grid=element_blank(),
             legend.position = "none"))

dev.off()



################################################################################




## let's filter out datasets with median of cells per patient <5 and number of patients <30

datasets_to_keep <- unlist(lapply(all_datasets, function(x){

  (median(c(x[[2]][1,]))>=5)&(ncol(x[[2]])>=30)
  
}))



filtered_datasets <- all_datasets[datasets_to_keep]

dataset_leaf <- sapply(names(datasets_to_keep), function(x){
  
  grepl("leaf",x)
  
})

sum(datasets_to_keep)

## leaves filtered out
nleaves(cluster_dendro)-sum(datasets_to_keep&dataset_leaf)

## nodes filtered out
nnodes(cluster_dendro)-nleaves(cluster_dendro)-sum(datasets_to_keep&(!dataset_leaf))

## lets save bed files

for (i in 1:length(filtered_datasets)){
  
  data_step <- filtered_datasets[[i]][[1]]
  
  nonzero_expression <- apply(data_step[,7:ncol(data_step)],1,
                              function(x){
                                
                                sum(x>0)/length(x)
                                
                              })
  names(nonzero_expression) <- data_step$phenotype_id
  nonzero_expression <- nonzero_expression[nonzero_expression>0.1]
  data_step <- data_step[data_step$phenotype_id%in%names(nonzero_expression),]
  
  mean_percent <- apply(data_step[,7:ncol(data_step)],2,
                        function(x){x/sum(x)})
  
  mean_percent <- rowMeans(mean_percent)
  names(mean_percent) <- data_step$phenotype_id
  mean_percent <- mean_percent[mean_percent>5e-7]
  
  data_step <- data_step[data_step$phenotype_id%in%names(mean_percent),]
  
  data_step[,1] <- paste0("chr",data_step[,1])
  
  write.table(data_step,
              paste0("bed/",names(filtered_datasets)[i],".bed"),
              quote = F, row.names = F, sep = "\t")
  
}

################################################################################
#########################_MAKE_METADATA_FILE_###################################

PCA_data <- read.delim(paste0("PCA.eigenvec"))
colnames(PCA_data) <- c("SampleID",paste0("genotype_",colnames(PCA_data)[2:6]))

meta <- read.delim("meta.tsv")
meta <- meta[meta$SampleID%in%PCA_data$SampleID,]
meta <- merge(meta,PCA_data, by = "SampleID")
meta$sex <- as.numeric(as.character(sapply(meta$sex, function(x){if (x=="M"){1}else{0}})))

write.csv(meta,"Covariates_biopsy.csv", quote = F, row.names = F)




