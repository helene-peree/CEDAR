ID_to_SNPname = function(ID){
  coordinates = data.frame(do.call("rbind", strsplit(ID, ":")))
  names(coordinates) = c("chr", "pos", "ref", "alt")
  
  SNPname = paste(coordinates[, "chr"], coordinates[, "pos"],
                  do.call(pmin, coordinates[, c("ref", "alt")]),
                  do.call(pmax, coordinates[, c("ref", "alt")]), sep = ":")
  
  return(SNPname)
}

read_ensembl = function(dataset){
  
  if (dataset == "blood"){
    file_name = paste(DocsDir, "ensembl_annotation_GRCh38_v105_start_end.tsv", sep = "/")
    
  } else if (dataset == "biopsies"){
    file_name = paste(DocsDir, "ensembl_annotation_GRCh38_v103_start_end.tsv", sep = "/")
  }
  
  ensembl_info = read.table(file_name, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ensembl_info[, "gene_name"] = ifelse(ensembl_info[, "external_gene_name"] != "", ensembl_info[, "external_gene_name"], ensembl_info[, "ensembl_gene_id"])
  
  return(ensembl_info)
}

read_genotypes = function(dataset){
  
  if (dataset == "blood"){
    file_name = paste(HDir, "raw_data", paste0("CEDAR-clean-snps-4-", chr, ".rds"), sep = "/")
    
  } else if (dataset == "biopsies"){
    file_name = paste(YDir, "raw_data", paste0("CEDAR_sc-", chr, ".rds"), sep = "/")
  }
  
  genotypes = readRDS(file_name)
  
  return(genotypes)
}

read_expression = function(dataset, cell_type, step){
  
  if (dataset == "blood"){
    RawDir = paste(HDir, "raw_data", sep = "/")
    ResultsDir = paste(HDir, "results", sep = "/")
    
    sample_ID = paste("ENSG", cell_type, "none", sep = "_")
    best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
    
  } else if (dataset == "biopsies") {
    RawDir = paste(YDir, "raw_data", sep = "/")
    ResultsDir = paste(YDir, "results", sep = "/")
    
    sample_ID = cell_type
    best = "full_peer"
  }
  
  if (step == "raw"){
    file_name = paste(RawDir, paste0(sample_ID, ".bed"), sep = "/")
    
  } else if (step == "normalized"){
    file_name = paste(ResultsDir, paste0("normalized_expression_", sample_ID),
                      paste0("expression_normalized_", sample_ID, ".bed"), sep = "/")
    
  } else if (step == "residualised"){
    file_name = paste(ResultsDir, paste0("residualised_expression_", sample_ID, "_", best),
                      paste0("expression_residualised_", sample_ID), sep = "/")
    
  } else if (step == "transformed"){
    file_name = paste(ResultsDir, paste0("residualised_expression_", sample_ID, "_", best),
                      paste0("expression_rank_normal_transformed_", sample_ID, ".bed"), sep = "/")
  }
  
  columns = read.table(file_name, comment.char = "", nrows = 1)
  
  expression = read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "",
                          colClasses = c("character", rep("NULL", 2), "character", rep("NULL", 2), rep("numeric", ncol(columns) - 6)))
  
  names(expression)[names(expression) %in% c("phenotype_id", "id")] = "gene"
  expression = expression[expression[, "#chr"] == chr,]
  expression = cbind(cell_type, expression)
  
  return(expression)
}

################################################################

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

if (analysis == "gene_specific"){
  iter = ""
  
} else if (analysis == "regulatory"){
  iter = paste0("_", iteration)
}

DocsDir = paste(prefix, "Analysis/RNAseq/docs", sep = "/")
DiseaseDir = paste(prefix, "Analysis/Theta_correlation/Pipeline_test/docs", sep = "/")
ViewerDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/PLATFORMS/GEN/Website", sep = "/")
ThetaDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta", dataset, sep = "/")
MetaDir = paste(ThetaDir, paste0("meta_analysis_", analysis, iter), sep = "/")
HDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA", sep = "/")
YDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap", sep = "/")

if (dataset == "blood"){
  eQTLDir = HDir
  
  hue_cols_8 = c("#d2a8ff", "#a371f7", "#6e40c9", "#cae8ff", "#a5d6ff", "#79c0ff", "#58a6ff", "#388bfd", "#1f6feb",
                 "#1158c7", "#0d419d", "#0c2d6b", "#051d4d", "#7ee787", "#3fb950", "#238636", "#d29922", "#9e6a03",
                 "#ffc680", "#f0883e", "#bd561d", "#ffa198", "#f85149", "#6e7681", "#ffbedd", "#f778ba", "#bf4b8a")
  
} else if (dataset == "biopsies"){
  eQTLDir = YDir
  hue_cols_8 = rep("black", 401)
  
} else if (dataset == "both"){
  eQTLDir = c(HDir, YDir)
}

RawDir = paste(eQTLDir, "raw_data", sep = "/")
ResultsDir = paste(eQTLDir, "results", sep = "/")

directories = list.files(path = ResultsDir, pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
directories = gsub("cis_eQTLs_permutational_|ENSG_|none_", "", directories)

overview = data.frame(cell_type = sub("_peer_.*", "", directories),
                      best_number_factors = sub(".*_peer_", "", directories))

cell_types = overview[, "cell_type"]

if (dataset == "blood"){
  cell_list = read.table(paste(DocsDir, "Cell_types_list.txt", sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  cell_list = cell_list[order(cell_list[, "new_order"]),]
  cell_type_names = cell_list[, "new_cell_type"]
  cell_type_old_names = cell_list[, "cell_type"]
  
} else if (dataset == "biopsies"){
  cell_type_names = sort(cell_types)
  cell_type_old_names = sort(cell_types)
}
