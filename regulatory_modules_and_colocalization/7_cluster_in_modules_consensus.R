library(qvalue)
library(ggplot2)

read_permutational = function(cell_type){
  
  if (! grepl("_pos_", cell_type)){
    best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
    
    file_name = paste(HDir, "results",
                      paste0("cis_eQTLs_permutational_ENSG_", cell_type, "_none_peer_", best),
                      paste0("cis.permutational_ENSG_", cell_type, "_none.txt.gz"), sep = "/")
    
  } else {
    file_name = paste(YDir, "results",
                      paste0("cis_eQTLs_permutational_", cell_type, "_peer_auto"),
                      paste0("cis.permutational_", cell_type, ".txt.gz"), sep = "/")
  }
  
  permutational = read.table(file_name, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
                             colClasses = c("character", rep("NULL", 4), "integer", "NULL", "character", "NULL", "integer", rep("NULL", 6), "numeric", rep("NULL", 3), "numeric"))
  
  names(permutational) = c("gene", "n_variants", "ID", "pos", "r_squared", "adj_beta_pval")
  permutational = permutational[permutational[, "n_variants"] >= 1,]
  permutational = cbind(cell_type, permutational)
}

annotate_modules = function(module){
  nr_EAPs = length(EAPs[EAPs[, "module"] == module, "EAP"])
  nr_genes = length(unique(EAPs[EAPs[, "module"] == module, "gene"]))
  cell_types = unique(EAPs[EAPs[, "module"] == module, "cell_type"])
  nr_cell_types = length(cell_types)
  vector = paste(ifelse(cell_types_old %in% cell_types, 1, 0), collapse = "")
  
  data.frame(module, nr_EAPs, nr_cell_types, nr_genes, vector)
}

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

dataset = "blood"

HDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA", sep = "/")
directories = list.files(path = paste(HDir, "results", sep = "/"), pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
directories = gsub("cis_eQTLs_permutational_ENSG_|none_", "", directories)

overview = data.frame(cell_type = sub("_peer_.*", "", directories),
                      best_number_factors = sub(".*_peer_", "", directories))

cell_types = overview[, "cell_type"]
cell_list = read.table(paste(prefix, "Analysis/RNAseq/docs/Cell_types_list.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cell_types_old = cell_list[order(cell_list[, "new_order"]), "cell_type"]

permutational = do.call("rbind", lapply(cell_types, read_permutational))
permutational[, "EAP"] = paste(permutational[, "gene"], permutational[, "cell_type"], sep = "_")

for (cell_type in cell_types){
  permutational[permutational[, "cell_type"] == cell_type, "qval"] = qvalue(permutational[permutational[, "cell_type"] == cell_type, "adj_beta_pval"])$qvalue
}

permutational[, "signif"] = ! is.na(permutational[, "qval"]) & permutational[, "qval"] <= 0.05

ThetaDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta", dataset, sep = "/")
EAPs_gene_specific = readRDS(paste(ThetaDir, paste0(dataset, "_EAPs_representative_gene_specific.rds"), sep = "/"))
EAPs_gene_specific[, "nr_cell_types"] = EAPs_gene_specific[, "nr_EAPs"]

#############################

if (! file.exists(paste(ThetaDir, "blood_EAPs_consensus.rds", sep = "/"))){
  thetas1 = readRDS(paste(ThetaDir, paste0(dataset, "_thetas_EAP_above_threshold.rds"), sep = "/"))
  thetas1 = thetas1[thetas1[, "gene_X"] == thetas1[, "gene_Y"],]
  thetas2 = readRDS(paste(ThetaDir, paste0(dataset, "_thetas_EAP_manual.rds"), sep = "/"))
  thetas2 = thetas2[thetas2[, "gene_X"] == thetas2[, "gene_Y"],]
  thetas3 = readRDS(paste(ThetaDir, paste0(dataset, "_thetas_EAP_consensus.rds"), sep = "/"))
  thetas3 = thetas3[thetas3[, "gene_X"] != thetas3[, "gene_Y"],]
  
  thetas = rbind(data.frame(origin = "thetas1", thetas1[, c("cell_type_X", "gene_X", "cell_type_Y", "gene_Y", "theta")]),
                 data.frame(origin = "thetas2", thetas2[, c("cell_type_X", "gene_X", "cell_type_Y", "gene_Y", "theta")]),
                 data.frame(origin = "thetas3", thetas3[, c("cell_type_X", "gene_X", "cell_type_Y", "gene_Y", "theta")]))
  
  thetas = thetas[abs(thetas[, "theta"]) >= 0.6,]
  thetas[, "EAP_X"] = paste(thetas[, "gene_X"], thetas[, "cell_type_X"], sep = "_")
  thetas[, "EAP_Y"] = paste(thetas[, "gene_Y"], thetas[, "cell_type_Y"], sep = "_")
  
  EAPs = data.frame(signif = TRUE, EAP = permutational[permutational[, "signif"], "EAP"], module = 0, stringsAsFactors = FALSE)
  
  for (EAP in EAPs[, "EAP"]){
    matched_EAPs = c(thetas[thetas[, "EAP_X"] == EAP, "EAP_Y"], thetas[thetas[, "EAP_Y"] == EAP, "EAP_X"])
    modules = unique(EAPs[EAPs[, "EAP"] %in% c(EAP, matched_EAPs), "module"])
  
    if (length(modules) == 1){
  
      if (modules == 0){
        EAPs[EAPs[, "EAP"] %in% c(EAP, matched_EAPs), "module"] = max(EAPs[, "module"])+ 1
      }
  
    } else if (length(modules) > 1){
      EAPs[EAPs[, "EAP"] %in% c(EAP, matched_EAPs), "module"] = min(setdiff(modules, 0))
  
      for (module in setdiff(modules, 0)){
        EAPs[EAPs[, "module"] == module, "module"] = min(setdiff(modules, 0))
      }
    }
  }
  
  saveRDS(EAPs, paste(ThetaDir, "blood_EAPs_consensus.rds", sep = "/"))
  
} else {
  EAPs = readRDS(paste(ThetaDir, "blood_EAPs_consensus.rds", sep = "/"))
}

EAPs[, "cell_type"] = sub(".*?_", "", EAPs[, "EAP"])
EAPs[, "gene"] = sub("_.*", "", EAPs[, "EAP"])

modules_freq = data.frame(table(EAPs[, "module"]))
names(modules_freq) = c("module", "nr_EAPs")
modules_freq = modules_freq[order(modules_freq[, "nr_EAPs"], decreasing = TRUE),]
modules_freq[, "new_module"] = 1:nrow(modules_freq)

EAPs = merge(EAPs, modules_freq[, c("module", "new_module", "nr_EAPs")], by = "module")
EAPs = EAPs[, setdiff(names(EAPs), "module")]
names(EAPs)[names(EAPs) == "new_module"] = "module"
rownames(EAPs) = NULL

modules = do.call("rbind", lapply(sort(unique(EAPs[, "module"])), annotate_modules))
saveRDS(modules, paste(ThetaDir, "blood_modules_consensus.rds", sep = "/"))
