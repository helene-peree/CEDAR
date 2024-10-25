args = commandArgs(trailingOnly = TRUE)
dataset = args[1] # "blood" or "biopsies" or "both"
analysis = args[2] # "gene_specific" or "regulatory"
iteration = args[3] # "first" or "second"
signif_only = FALSE #TRUE or FALSE

################################################################

library(qvalue)

################################################################

source("function_definitions.R")

read_ensembl = function(dataset){
  
  if (dataset == "blood"){
    file_name = paste(prefix, "Analysis/RNAseq/docs/ensembl_annotation_GRCh38_v105_start_end.tsv", sep = "/")
    
  } else if (dataset == "biopsies") {
    file_name = paste(prefix, "Analysis/RNAseq/docs/ensembl_annotation_GRCh38_v103_start_end.tsv", sep = "/")
  }
  
  ensembl_info = read.table(file_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ensembl_info[, "gene_name"] = ifelse(ensembl_info[, "external_gene_name"] != "", ensembl_info[, "external_gene_name"], ensembl_info[, "ensembl_gene_id"])
  ensembl_info[, "dataset_gene"] = paste(dataset, ensembl_info[, "ensembl_gene_id"], sep = "_")
  ensembl_info = ensembl_info[, c("dataset_gene", "ensembl_gene_id", "gene_name", "chromosome_name", "transcription_start_site", "start_position", "end_position", "strand")]
  
  return(ensembl_info)
}

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

##########

threshold_blood = 0.9
threshold_biopsies = 1.9

if (analysis == "gene_specific"){
  iter = ""
  
} else if (analysis == "regulatory"){
  iter = paste0("_", iteration)
}

DocsDir = paste(prefix, "Analysis/RNAseq/docs", sep = "/")
ScriptsDir = paste(prefix, "Analysis/eQTL/scripts", sep = "/")
HclustDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/hclust", sep = "/")
ThetaDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta", sep = "/")
HDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA", sep = "/")
YDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap", sep = "/")

EAPs_wo_filename = paste(ThetaDir, dataset, paste0(dataset, "_EAPs_", analysis, iter, ".rds"), sep = "/")
EAPs_wi_filename = paste(ThetaDir, dataset, paste0(dataset, "_EAPs_representative_", analysis, iter, ".rds"), sep = "/")

if (signif_only){
  modules_filename = paste(ThetaDir, dataset, paste0(dataset, "_modules_", analysis, iter, "_signif.rds"), sep = "/")
  
} else {
  modules_filename = paste(ThetaDir, dataset, paste0(dataset, "_modules_", analysis, iter, ".rds"), sep = "/")
}

if (dataset != "both"){
  ensembl_info = read_ensembl(dataset)
  
} else {
  ensembl_info = rbind(read_ensembl("blood"),
                       read_ensembl("biopsies"))
}

directories = list.files(path = paste(HDir, "results", sep = "/"), pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
directories = gsub("cis_eQTLs_permutational_ENSG_|none_", "", directories)

overview = data.frame(cell_type = sub("_peer_.*", "", directories),
                      best_number_factors = sub(".*_peer_", "", directories))

cell_types_blood = overview[, "cell_type"]
cell_list = read.table(paste(DocsDir, "Cell_types_list.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
cell_types_blood_old = cell_list[order(cell_list[, "new_order"]), "cell_type"]
cell_types_blood_new = cell_list[order(cell_list[, "new_order"]), "new_cell_type"]

directories = list.files(path = paste(YDir, "results", sep = "/"), pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
cell_types_biopsies = gsub("cis_eQTLs_permutational_|_peer_auto", "", directories)

if (dataset == "blood"){
  cell_types = cell_types_blood
  cell_types_old = cell_types_blood_old
  cell_types_new = cell_types_blood_new
  
} else if (dataset == "biopsies"){
  cell_types = cell_types_biopsies
  cell_types_old = cell_types_biopsies
  cell_types_new = cell_types_biopsies
  
} else if (dataset == "both"){
  cell_types = c(cell_types_blood, cell_types_biopsies)
  cell_types_old = c(cell_types_blood_old, cell_types_biopsies)
  cell_types_new = c(cell_types_blood_new, cell_types_biopsies)
}

permutational = do.call("rbind", lapply(cell_types, read_permutational))
permutational[, "EAP"] = paste(permutational[, "gene"], permutational[, "cell_type"], sep = "_")

for (cell_type in cell_types){
  permutational[permutational[, "cell_type"] == cell_type, "qval"] = qvalue(permutational[permutational[, "cell_type"] == cell_type, "adj_beta_pval"])$qvalue
}

permutational[, "signif"] = ! is.na(permutational[, "qval"]) & permutational[, "qval"] <= 0.05

########## Import thetas

thetas_files_blood = paste(ThetaDir, "blood", c("blood_thetas_EAP.rds", "blood_thetas_EAP_manual.rds"), sep = "/")
thetas_files_biopsies = paste(ThetaDir, "biopsies", c("biopsies_thetas_EAP_above_threshold.rds", "biopsies_thetas_EAP_manual.rds"), sep = "/")
thetas_files_both = paste(ThetaDir, "both", c("both_thetas_EAP_above_threshold.rds", "both_thetas_EAP_manual.rds"), sep = "/")

if (dataset == "blood"){
  thetas_files = thetas_files_blood
  
} else if (dataset == "biopsies"){
  thetas_files = thetas_files_biopsies
  
} else if (dataset == "both"){
  thetas_files = c(thetas_files_blood, thetas_files_biopsies, thetas_files_both)
}

thetas = do.call("rbind", lapply(thetas_files, function(x){
  x = readRDS(x)
  x = x[abs(x[, "theta"]) >= 0.6, c("cell_type_X", "gene_X", "cell_type_Y", "gene_Y", "theta", "n_variants")]
}))

thetas[, "EAP_X"] = paste(thetas[, "gene_X"], thetas[, "cell_type_X"], sep = "_")
thetas[, "EAP_Y"] = paste(thetas[, "gene_Y"], thetas[, "cell_type_Y"], sep = "_")
thetas = thetas[!is.na(thetas[, "theta"]),]

if (analysis == "gene_specific"){
  thetas = thetas[thetas[, "gene_X"] == thetas[, "gene_Y"],]
  
} else if (analysis == "regulatory" & iteration = "second"){
  OutputDir = paste(ThetaDir, dataset, paste0("check_", analysis, "_first"), sep = "/")
  files = list.files(path = OutputDir, pattern = "theta_chr", full.names = TRUE)
  coloc = do.call("rbind", lapply(files, read.table, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
  coloc[, "EAP"] = paste(coloc[, "gene"], coloc[, "cell_type"], sep = "_")
  kept = coloc[abs(coloc[, "theta"]) >= 0.6, "EAP"]
  not_kept = coloc[abs(coloc[, "theta"]) < 0.6, "EAP"]
  
  thetas[thetas[, "EAP_X"] %in% kept & thetas[, "EAP_Y"] %in% not_kept, "theta"] = 0
  thetas[thetas[, "EAP_X"] %in% not_kept & thetas[, "EAP_Y"] %in% kept, "theta"] = 0
  thetas = thetas[abs(thetas[, "theta"]) >= 0.6,]
}

########## Assign each EAP to a module

EAPs_significant = data.frame(EAP = permutational[permutational[, "signif"], "EAP"], module = 0, stringsAsFactors = FALSE)
EAPs_non_significant = data.frame(EAP = permutational[! permutational[, "signif"] & ! permutational[, "signif"] & ((! grepl("_pos_", permutational[, "cell_type"]) & -log10(permutational[, "adj_beta_pval"]) >= threshold_blood) | (grepl("_pos_", permutational[, "cell_type"]) & -log10(permutational[, "adj_beta_pval"]) >= threshold_biopsies)) & permutational[, "EAP"] %in% unique(c(thetas[, "EAP_X"], thetas[, "EAP_Y"])), "EAP"], module = 0, stringsAsFactors = FALSE)

for (EAP in EAPs_significant[, "EAP"]){
  matched_EAPs = c(thetas[thetas[, "EAP_X"] == EAP, "EAP_Y"], thetas[thetas[, "EAP_Y"] == EAP, "EAP_X"])
  modules = unique(EAPs_significant[EAPs_significant[, "EAP"] %in% c(EAP, matched_EAPs), "module"])
  
  if (length(modules) == 1){
    
    if (modules == 0){
      EAPs_significant[EAPs_significant[, "EAP"] %in% c(EAP, matched_EAPs), "module"] = max(EAPs_significant[, "module"])+ 1
    }
    
  } else if (length(modules) > 1){
    EAPs_significant[EAPs_significant[, "EAP"] %in% c(EAP, matched_EAPs), "module"] = min(setdiff(modules, 0))
    
    for (module in setdiff(modules, 0)){
      EAPs_significant[EAPs_significant[, "module"] == module, "module"] = min(setdiff(modules, 0))
    }
  }
}

for (EAP in EAPs_non_significant[, "EAP"]){
  matched_EAPs = thetas[thetas[, "EAP_X"] == EAP | thetas[, "EAP_Y"] == EAP,]
  matched_EAPs = matched_EAPs[matched_EAPs[, "EAP_X"] %in% permutational[permutational[, "signif"], "EAP"] | matched_EAPs[, "EAP_Y"] %in% permutational[permutational[, "signif"], "EAP"],]
  
  if (nrow(matched_EAPs) > 0){
    matched_EAPs = matched_EAPs[which.max(abs(matched_EAPs[, "theta"])), c("EAP_X", "EAP_Y")]
    module = EAPs_significant[EAPs_significant[, "EAP"] %in% matched_EAPs, "module"]
    EAPs_non_significant[EAPs_non_significant[, "EAP"] == EAP, "module"] = module
    
  } else{
    EAPs_non_significant = EAPs_non_significant[EAPs_non_significant[, "EAP"] != EAP,]
  }
}

EAPs = rbind(data.frame(signif = TRUE, EAPs_significant, stringsAsFactors = FALSE),
             data.frame(signif = FALSE, EAPs_non_significant, stringsAsFactors = FALSE))

EAPs[, "cell_type"] = sub(".*?_", "", EAPs[, "EAP"])
EAPs[, "gene"] = sub("_.*", "", EAPs[, "EAP"])
EAPs[, "dataset_gene"] = ifelse(! grepl("_pos_", EAPs[, "cell_type"]), paste0("blood_", EAPs[, "gene"]), paste0("biopsies_", EAPs[, "gene"]))

modules_freq = data.frame(table(EAPs[, "module"]))
names(modules_freq) = c("module", "nr_EAPs")
modules_freq = modules_freq[order(modules_freq[, "nr_EAPs"], decreasing = TRUE),]
modules_freq[, "new_module"] = 1:nrow(modules_freq)

EAPs = merge(EAPs, modules_freq[, c("module", "new_module", "nr_EAPs")], by = "module")
EAPs = merge(EAPs, permutational[, c("EAP", "ID", "pos", "adj_beta_pval")], by = "EAP")
EAPs = merge(EAPs, ensembl_info[, c("dataset_gene", "ensembl_gene_id", "gene_name", "chromosome_name", "transcription_start_site", "start_position", "end_position", "strand")], by = "dataset_gene")
EAPs = EAPs[order(EAPs[, "transcription_start_site"]),]
EAPs = EAPs[order(EAPs[, "chromosome_name"]),]
EAPs = EAPs[, setdiff(names(EAPs), "module")]
names(EAPs)[names(EAPs) == "new_module"] = "module"
rownames(EAPs) = NULL

saveRDS(EAPs, EAPs_wo_filename)

########## Summary for each module

annotate_modules = function(module){
  
  if (signif_only){
    EAPs_subset = EAPs[EAPs[, "module"] == module & EAPs[, "signif"],]
    
  } else {
    EAPs_subset = EAPs[EAPs[, "module"] == module,]
  }
  
  thetas_subset = thetas[thetas[, "EAP_X"] %in% EAPs_subset[, "EAP"] & thetas[, "EAP_Y"] %in% EAPs_subset[, "EAP"],]
  
  EAPs_signif_subset = EAPs_subset[EAPs_subset[, "signif"],]
  representative = EAPs_signif_subset[which.min(EAPs_signif_subset[, "adj_beta_pval"]), "EAP"]
  nr_EAPs = length(EAPs_subset[, "EAP"])
  nr_signif_EAPs = length(EAPs_subset[EAPs_subset[, "EAP"] %in% permutational[permutational[, "signif"], "EAP"], "EAP"])
  nr_non_signif_EAPs = length(EAPs_subset[EAPs_subset[, "EAP"] %in% permutational[! permutational[, "signif"], "EAP"], "EAP"])
  cell_types = unique(EAPs_subset[, "cell_type"])
  blood_cell_types = cell_types[! grepl("_pos_", cell_types)]
  gut_nodes_or_leaves = cell_types[grepl("_pos_", cell_types)]
  nr_cell_types = length(cell_types)
  nr_blood_cell_types = length(blood_cell_types)
  nr_gut_nodes_or_leaves = length(gut_nodes_or_leaves)
  nr_datasets = sum(any(grepl("_pos_", cell_types))) + sum(any(! grepl("_pos_", cell_types)))
  nr_locations = length(unique(sapply(grep("leaf", EAPs_subset[, "EAP"], value = TRUE), function(x){strsplit(as.character(x), "_")[[1]][3]})))
  
  vector = paste(ifelse(cell_types_old %in% cell_types, 1, 0), collapse = "")
  
  genes = unique(EAPs_subset[, "gene"])
  nr_genes = length(genes)
  gene_ids = paste(sort(genes), collapse = ", ")
  gene_names = paste(sort(unique(EAPs_subset[, "gene_name"])), collapse = ", ")
  
  nr_links = nrow(thetas_subset)
  nr_pos_links = nrow(thetas_subset[sign(thetas_subset[, "theta"]) == 1,])
  nr_neg_links = nrow(thetas_subset[sign(thetas_subset[, "theta"]) == -1,])
  max_nr_links = (nr_EAPs * (nr_EAPs - 1))/2
  
  if (nr_EAPs > 1){
    proportion_links = round(nr_links / max_nr_links, 2)
    mean_abs_theta = round(mean(abs(thetas_subset[, "theta"])), 2)
    
  } else {
    proportion_links = NA
    mean_abs_theta = NA
  }
  
  if (dataset != "both"){
    left_border = min(EAPs_subset[, "transcription_start_site"]) - 10^6
    right_border = max(EAPs_subset[, "transcription_start_site"]) + 10^6
    chromosome_name = unique(EAPs_subset[, "chromosome_name"])
    
  } else {
    left_border = 0
    right_border = 0
    chromosome_name = 0
  }
  
  data.frame(module, nr_EAPs, nr_signif_EAPs, nr_non_signif_EAPs, nr_links,
             nr_pos_links, nr_neg_links, max_nr_links, proportion_links,
             mean_abs_theta, nr_cell_types, nr_blood_cell_types, nr_gut_nodes_or_leaves,
             nr_datasets, nr_locations, nr_genes, gene_ids, gene_names,
             representative, chromosome_name, left_border, right_border, vector)
}

modules = do.call("rbind", lapply(sort(unique(EAPs[, "module"])), annotate_modules))
rownames(modules) = NULL

for (cell_type in cell_types_new){
  modules[, cell_type] = sapply(modules[, "vector"], function(x){strsplit(as.character(x), "")[[1]][cell_types_new == cell_type]})
}

saveRDS(modules, modules_filename)

########## Assign up and down signs to each EAP

EAPs[, "known"] = ifelse(EAPs[, "EAP"] %in% modules[, "representative"], "representative", NA)
EAPs[, "theta_sign"] = ifelse(EAPs[, "EAP"] %in% modules[, "representative"], 1, NA)
EAPs[, "i"] = ifelse(EAPs[, "EAP"] %in% modules[, "representative"], 0, NA)
i = 0

while(any(is.na(EAPs[, "theta_sign"]))){
  i = i + 1
  
  for (EAP in EAPs[is.na(EAPs[, "theta_sign"]), "EAP"]){
    module = EAPs[EAPs[, "EAP"] == EAP, "module"]
    known_EAPs = EAPs[EAPs[, "module"] == module & EAPs[, "i"] == i - 1 & !is.na(EAPs[, "theta_sign"]), "EAP"]
    thetas_subset = thetas[thetas[, "EAP_X"] %in% c(EAP, known_EAPs) & thetas[, "EAP_Y"] %in% c(EAP, known_EAPs),]
    thetas_subset = thetas_subset[thetas_subset[, "EAP_X"] == EAP | thetas_subset[, "EAP_Y"] == EAP,]
    thetas_subset = thetas_subset[thetas_subset[, "EAP_X"] %in% permutational[permutational[, "signif"], "EAP"] | thetas_subset[, "EAP_Y"] %in% permutational[permutational[, "signif"], "EAP"],]
    thetas_subset = thetas_subset[which.max(abs(thetas_subset[, "theta"])),]
    
    if (nrow(thetas_subset) > 0){
      known = as.character(setdiff(thetas_subset[, c("EAP_X", "EAP_Y")], EAP))
      previous_sign = EAPs[EAPs[, "EAP"] == known, "theta_sign"]
      current_sign = sign(thetas_subset[, "theta"])
      EAPs[EAPs[, "EAP"] == EAP, "known"] = known
      EAPs[EAPs[, "EAP"] == EAP, "theta_sign"] = previous_sign * current_sign
      EAPs[EAPs[, "EAP"] == EAP, "i"] = i
    }
  }
}

saveRDS(EAPs, EAPs_wi_filename)
