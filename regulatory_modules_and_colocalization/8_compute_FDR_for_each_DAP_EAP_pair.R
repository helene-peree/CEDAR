args = commandArgs(trailingOnly = TRUE)
dataset = args[1] # "blood" or "biopsies"

################################################################

library(qvalue)

################################################################

read_permutational = function(cell_type, overview, eQTLDir){
  
  if (! grepl("_pos_", cell_type)){
    sample_ID = paste("ENSG", cell_type, "none", sep = "_")
    
  } else {
    sample_ID = cell_type
  }
  
  best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
  
  file_name = paste(eQTLDir,
                    paste0("cis_eQTLs_permutational_", sample_ID, "_peer_", best),
                    paste0("cis.permutational_", sample_ID, ".txt.gz"), sep = "/")
  
  permutational = read.table(file_name, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  
  names(permutational) = c("phe_id", "phe_chr", "phe_from", "phe_to", "phe_strd",
                           "n_var_in_cis", "dist_phe_var", "var_id", "var_chr", "var_from",
                           "var_to", "dof1", "dof2", "bml1", "bml2",
                           "nom_pval", "r_squared", "slope", "slope_se", "adj_emp_pval",
                           "adj_beta_pval")
  
  permutational = cbind(cell_type, permutational)
}

prepare_tables = function(ThetaDir, eQTLDir){
  files = list.files(path = ThetaDir, pattern = ".tsv", full.names = TRUE, recursive = TRUE)
  coloc = do.call("rbind", lapply(files, function(file_name){x = read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)}))
  coloc[coloc[, "p_value"] == 0, "p_value"] = 1 / 10000
  directories = list.files(path = eQTLDir, pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
  directories = gsub("cis_eQTLs_permutational_|ENSG_|none_", "", directories)
  overview = data.frame(cell_type = sub("_peer_.*", "", directories), best_number_factors = sub(".*_peer_", "", directories))
  cell_types = overview[, "cell_type"]
  permutational = do.call("rbind", lapply(cell_types, read_permutational, overview, eQTLDir))
  
  for (cell_type in cell_types){
    permutational[permutational[, "cell_type"] == cell_type, "qval"] = qvalue(permutational[permutational[, "cell_type"] == cell_type, "adj_beta_pval"])$qvalue
  }
  
  permutational[, "signif"] = ! is.na(permutational[, "qval"]) & permutational[, "qval"] <= 0.05
  permutational[, "eGene"] = permutational[, "phe_id"] %in% permutational[permutational[, "signif"], "phe_id"]
  coloc = merge(coloc, permutational[, c("cell_type", "phe_id", "nom_pval", "adj_beta_pval", "qval", "signif", "eGene")], by.x = c("cell_type", "gene"), by.y = c("cell_type", "phe_id"))
}

################################################################

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

if (dataset == "blood"){
  real = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/blood/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                        paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA/results", sep = "/"))
  
  permu = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/blood_permu/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                         paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_permu/results", sep = "/"))
  
} else if (dataset == "biopsies"){
  real = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                        paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap/results", sep = "/"))
  
  permu_1 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_1/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_1/results", sep = "/"))
  
  permu_2 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_2/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_2/results", sep = "/"))
  
  permu_3 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_3/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_3/results", sep = "/"))
  
  permu_4 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_4/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_4/results", sep = "/"))
  
  permu_5 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_5/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_5/results", sep = "/"))
  
  permu_6 = prepare_tables(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/biopsies_permu_6/DAP/10000_permutations_0.1_manual_0", sep = "/"),
                                    paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/CEDAR_sc_full_shuffle/CEDAR_sc_qtl_6/results", sep = "/"))
  
  permu = rbind(permu_1, permu_2, permu_3, permu_4, permu_5, permu_6)
  permu = permu[permu[, "cell_type"] %in% real[, "cell_type"],]
  rm(permu_1, permu_2, permu_3, permu_4, permu_5, permu_6)
}

real[, "abs_theta"] = abs(real[, "theta"])
real[, "log_theta_pval"] = -log10(real[, "p_value"])
real[, "log_eQTL_pval"] = -log10(real[, "adj_beta_pval"])
real_subset = real[real[, "n_variants"] >= 8 & real[, "abs_theta"] >= 0.6,]

permu[, "abs_theta"] = abs(permu[, "theta"])
permu[, "log_theta_pval"] = -log10(permu[, "p_value"])
permu[, "log_eQTL_pval"] = -log10(permu[, "adj_beta_pval"])
permu_subset = permu[permu[, "n_variants"] >= 8 & permu[, "abs_theta"] >= 0.6,]

for (i in 1:nrow(real)){
  theta = real[i, "abs_theta"]
  log10_theta_pval = real[i, "log_theta_pval"]
  log10_eQTL_pval = real[i, "log_eQTL_pval"]
  
  nb_pos_real_data = nrow(real_subset[real_subset[, "abs_theta"] >= theta & real_subset[, "log_theta_pval"] >= log10_theta_pval & real_subset[, "log_eQTL_pval"] >= log10_eQTL_pval,])
  nb_pos_permu_data = nrow(permu_subset[permu_subset[, "abs_theta"] >= theta & permu_subset[, "log_theta_pval"] >= log10_theta_pval & permu_subset[, "log_eQTL_pval"] >= log10_eQTL_pval,])
  
  real[i, "FDR"] = nb_pos_permu_data / (nb_pos_permu_data + nb_pos_real_data)
}

real[, "tier"] = ifelse(real[, "n_variants"] >= 8 & real[, "abs_theta"] >= 0.6 & real[, "FDR"] <= 0.05, "tier_1", ifelse(real[, "n_variants"] >= 8 & real[, "abs_theta"] >= 0.6 & real[, "FDR"] > 0.05 & real[, "FDR"] <= 0.1, "tier_2", "other"))
saveRDS(real, paste0("coloc_", dataset, "_FDR_threshold_all.rds"))
