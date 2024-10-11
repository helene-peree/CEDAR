library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(DescTools)
library(scales)
library(harmony)


path_scdata <- ""
samples_info <- read.csv("samples_info_rerun.csv")

samples_info <- samples_info[!samples_info$Sample%in%c("B2","B8","B15"),]

cell_statistics <- as.data.frame(matrix(ncol = 16, nrow = nrow(samples_info)))
colnames(cell_statistics) <- c("sample","n_Cells","mto_percent","mto_percent_filt","n_cells_filt",
                               "Doublet_percent_HTO_km","negative_percent_HTO_km",
                               "Doublet_percent_MULTI","negative_percent_MULTI",
                               "doubled_singlet","doublets_for_sure","negative_doublet",
                               "negative_for_sure","negative_singlet","singlets_same_type",
                               "singlets_wrong_type")

cell_labels <- vector(mode = "list", length = nrow(samples_info))
names(cell_labels) <- samples_info[1:nrow(samples_info),1]

seurat_objects <- vector(mode = "list", length = nrow(samples_info))
names(seurat_objects) <- samples_info[1:nrow(samples_info),1]

for (i in 1:nrow(samples_info)){

  print(samples_info[i,1])

  #### LOAD 10x data
  
  
  ##############################################################################

  
  feature_names <- read.delim(paste0(path_scdata,samples_info[i,2],
                                        "features.tsv.gz"), header = F)
  feature_names$V1 <- mapply(function(x,y){
    
    if (grepl("Hashtag",x)){x <- y} else {x <- paste0(x,"_",y)}
    x
  }, x = feature_names$V1, y = feature_names$V2)
  
  system(paste0("mv ", path_scdata,samples_info[i,2],
                "features.tsv.gz ",path_scdata,samples_info[i,2], 
                "features_old.tsv.gz"))
  write.table(feature_names,paste0(path_scdata,samples_info[i,2],
                                   "features.tsv"), col.names = F, row.names = F)
  system(paste0("gzip ", path_scdata,samples_info[i,2],"features.tsv"))
  

  ##############################################################################
  
  
  data_10x <- Read10X(paste0(path_scdata,samples_info[i,2]), gene.column = 1)
  
  colnames(data_10x$`Gene Expression`) <- 
    colnames(data_10x$`Antibody Capture`) <- paste0(colnames(data_10x$`Antibody Capture`),
                                                    "-",samples_info[i,1])

  if (samples_info[i,1]=="B33"|
      samples_info[i,1]=="B53"){rownames(data_10x$`Antibody Capture`)[10] <- "TC_LP.2"} else
        if (samples_info[i,1]=="B53"){
          rownames(data_10x$`Antibody Capture`)[c(5,6,8:10)] <- c("MIX_RE_IL_LP",
                                                                  "MIX_RE_IL_LP.1",
                                                                  "MIX_RE_IL_LP.2",
                                                                  "MIX_RE_IL_LP.3",
                                                                  "TC_LP.2")

        } else if (samples_info[i,1]%in%c("B26","B28","B39","B41","B44","B60")){

          data_10x$`Antibody Capture` <- data_10x$`Antibody Capture`[!(grepl("IL",rownames(data_10x$`Antibody Capture`))),]

        } else if (samples_info[i,1]%in%paste0("B", 2:10)) {

          rownames(data_10x$`Antibody Capture`)[c(5,6)] <- c("RE_EC","RE_IEL")

        } else {}


  #### Save the abundance of each barcode and No of barcodes counts per each cell

  barcode_numbers <- rowSums(data_10x$`Antibody Capture`)
  barcodes_cells_number <- colSums(data_10x$`Antibody Capture`)

  #### If any barcode has less than or equal to 300 of reads - remove it

  if (any(barcode_numbers<=300)){

    data_10x$`Antibody Capture` <- data_10x$`Antibody Capture`[barcode_numbers>300,]

  }

  cell_statistics$sample[i] <- samples_info[i,1]

  #### Create seurat object with gene expression data
  
  seurat_obj <- CreateSeuratObject(counts = data_10x$`Gene Expression`, project = samples_info[i,1])
  cell_statistics$n_Cells[i] <- length(seurat_obj$nCount_RNA)

  #### Add barcode counts to the object

  seurat_obj[['HTO']] <- CreateAssayObject(counts = data_10x$`Antibody Capture`)

  #### Compute % of mitochondrial reads per each cell

  seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-", col.name = "percent.mt")
  cell_statistics$mto_percent[i] <- mean(seurat_obj$percent.mt)

  #### Filter cells with % of mitochondrial reads more than or equal to 50
  #### and number of genes expressed less than 200

  seurat_obj <- subset(seurat_obj, subset = percent.mt < 50)
  seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200)

  #### Filter cells with less than 5 barcode counts

  seurat_obj <- seurat_obj[,!colnames(seurat_obj)%in%names(barcodes_cells_number[barcodes_cells_number<=5])]

  cell_statistics$mto_percent_filt[i] <- mean(seurat_obj$percent.mt)
  cell_statistics$n_cells_filt[i] <- length(seurat_obj$nCount_RNA)

  ##############_RUN_HTO_KMEANS_DEMULTIPLEXING_#################################

  #### make CLR normalization along features

  seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "CLR")

  #### Run demultiplexing - normalization is the same

  seurat_obj <- HTODemux(seurat_obj, assay = "HTO", positive.quantile = 0.999,
                         kfunc = "kmeans")

  #### Save cell classification results
  class_km <- as.factor(seurat_obj$HTO_classification.global)
  hash.ID_km <- as.character(seurat_obj$hash.ID)

  cell_statistics$Doublet_percent_HTO_km[i] <- c((table(seurat_obj$HTO_classification.global)/
                                                    sum(table(seurat_obj$HTO_classification.global)))*100)[1]
  cell_statistics$negative_percent_HTO_km[i] <- c((table(seurat_obj$HTO_classification.global)/
                                                     sum(table(seurat_obj$HTO_classification.global)))*100)[2]

  ##############################################################################

  class_km_doublets <- seurat_obj$HTO_classification

  class_km_doublets_mod <- sapply(class_km_doublets, function(x){

    gsub(".1","", fixed = T, gsub(".2","", fixed = T, x))

  })

  class_km_doublets_mod <- sapply(class_km_doublets_mod, function(x){

    y <- unlist(strsplit(x, "_", fixed = T))

    if (length(y)==1){y} else if (length(unique(y)) == 1){y[1]} else{paste(y,collapse="_")}

  })

  class_km_doublets_new <- sapply(class_km_doublets_mod, function(x){

    if(grepl("_",x)){"Doublet"}else{x}

  })

  ##############################################################################


  ###############_RUN_HTO_MULTI_DEMULTIPLEXING_#################################

  #### Normalize the data - log normalization along features

  seurat_obj <- NormalizeData(seurat_obj, assay = "HTO", normalization.method = "LogNormalize")

  #### Scale data

  seurat_obj <- ScaleData(seurat_obj, assay = "HTO")

  #### Make demultiplexing

  seurat_obj <- MULTIseqDemux(seurat_obj, assay = "HTO", autoThresh = T)

  #### Make same type of "global" cell classification for multi
  hash.ID_MULTI <- as.character(seurat_obj$MULTI_ID)

  class_multi <- seurat_obj$MULTI_ID
  levels(class_multi) <- sapply(levels(class_multi), function(x){

    if (x!="Doublet"&x!="Negative"){x <- "Singlet"} else {x}

  })

  if (!"Negative"%in%names(table(class_multi))){

    levels(class_multi) <- c("Doublet","Singlet","Negative")

  }

  class_multi <- factor(class_multi, levels = c("Doublet","Negative","Singlet"))

  cell_statistics$Doublet_percent_MULTI[i] <- c((table(class_multi)/
                                                   sum(table(class_multi)))*100)[1]
  cell_statistics$negative_percent_MULTI[i] <- c((table(class_multi)/
                                                    sum(table(class_multi)))*100)[2]

  ##############################################################################

  class_MULTI_doublets <- seurat_obj$MULTI_classification

  class_MULTI_doublets_mod <- sapply(class_MULTI_doublets, function(x){

    gsub(".1","", fixed = T, gsub(".2","", fixed = T, x))

  })

  class_MULTI_doublets_mod <- sapply(class_MULTI_doublets_mod, function(x){

    y <- unlist(strsplit(x, "_", fixed = T))

    if (length(y)==1){y} else if (length(unique(y)) == 1){y[1]} else{paste(y,collapse="_")}

  })

  class_MULTI_doublets_new <- sapply(class_MULTI_doublets_mod, function(x){

    if(grepl("_",x)|is.na(x)){"Doublet"}else{x}

  })

  ##############################################################################

  cells_for_sure <- mapply(function(x,y){

    if((x==y)&(x!="Doublet"&x!="Negative")&(y!="Doublet"&y!="Negative")){"singlets_same_type"} else
      if ((x!=y)&(x!="Doublet"&x!="Negative")&(y!="Doublet"&y!="Negative")){"singlets_wrong_type"} else
        if ((x=="Negative")&(y=="Negative")){"negative_for_sure"} else
          if ((x=="Doublet")&(y=="Doublet")){"doublets_for_sure"} else
            if (((x=="Negative")&(y!="Doublet"))|((y=="Negative")&(x!="Doublet"))){"negative_singlet"} else
              if (((x=="Doublet")&(y!="Negative"))|((y=="Doublet")&(x!="Negative"))){"doubled_singlet"} else
              {"negative_doublet"}

  }, x = class_km_doublets_new, y = class_MULTI_doublets_new)

  if (any(!colnames(cell_statistics[10:16])%in%names(c(table(cells_for_sure))))){

    cells_for_sure_res <- c(table(cells_for_sure))
    cells_for_sure_res <- cells_for_sure_res[match(colnames(cell_statistics[10:16]),
                                                   names(cells_for_sure_res))]
    cells_for_sure_res[is.na(cells_for_sure_res)] <- 0
    cell_statistics[i, 10:16] <- cells_for_sure_res

  } else {

    cell_statistics[i, 10:16] <- c(table(cells_for_sure))

  }

  seurat_obj@meta.data$cells_for_sure <- cells_for_sure

  cell_labels[[i]] <- cells_for_sure
  seurat_objects[[i]] <- seurat_obj

}

rm(seurat_obj, hash.ID_km, hash.ID_MULTI, class_km, class_km_doublets,
   class_km_doublets_mod, class_km_doublets_new, class_multi, class_MULTI_doublets,
   class_MULTI_doublets_mod, class_MULTI_doublets_new); gc()

################################################################################

saveRDS(cell_statistics,"cell_statistics.RDS")
saveRDS(seurat_objects, "seurat_objects_all.RDS")
saveRDS(cell_labels, "cell_labels_all.RDS")

