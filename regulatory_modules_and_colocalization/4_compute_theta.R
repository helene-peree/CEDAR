suppressMessages(library(Rfast))

################################################################

ID_to_SNPname = function(ID){
  coordinates = data.frame(do.call("rbind", strsplit(ID, ":")))
  names(coordinates) = c("chr", "pos", "ref", "alt")
  
  SNPname = paste(coordinates[, "chr"], coordinates[, "pos"],
                  do.call(pmin, coordinates[, c("ref", "alt")]),
                  do.call(pmax, coordinates[, c("ref", "alt")]), sep = ":")
  
  return(SNPname)
}

read_genotypes = function(cell_type, individuals){
  
  if (! grepl("pos", cell_type)){
    file_name = paste(WorkingDir, "raw_data", paste0("CEDAR-clean-snps-4-", chr, ".rds"), sep = "/")
    
  } else {
    file_name = paste(WorkingDir, "raw_data", paste0("CEDAR_sc-", chr, ".rds"), sep = "/")
  }
  
  genotypes = readRDS(file_name)
  genotypes[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", row.names(genotypes))))
  genotypes[, "MAF"] = rowSums(round(genotypes[, individuals])) / (2 * length(individuals))
  
  return(genotypes)
}

read_expression = function(cell_type){
  
  if (! grepl("pos", cell_type)){
    best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
    
    file_name = paste(WorkingDir, "results",
                      paste0("residualised_expression_ENSG_", cell_type, "_none_", best),
                      paste0("expression_rank_normal_transformed_ENSG_", cell_type, "_none.bed"), sep = "/")
    
  } else {
    file_name = paste(WorkingDir, "results",
                      paste0("residualised_expression_", cell_type, "_full_peer"),
                      paste0("expression_rank_normal_transformed_", cell_type, ".bed"), sep = "/")
  }
  
  columns = read.table(file_name, comment.char = "", nrows = 1)
  
  expression = read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "",
                          colClasses = c("character", rep("NULL", 2), "character", rep("NULL", 2), rep("numeric", ncol(columns) - 6)))
  
  names(expression)[names(expression) == "id"] = "gene"
  expression = expression[expression[, "#chr"] == chr,]
  expression = cbind(cell_type, expression)
  
  return(expression)
}

read_significant = function(cell_type){
  
  if (! grepl("pos", cell_type)){
    best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
    
    file_name = paste(WorkingDir, "results",
                      paste0("cis_eQTLs_permutational_ENSG_", cell_type, "_none_peer_", best),
                      paste0("cis.permutational_ENSG_", cell_type, "_none_FDR.significant.txt"), sep = "/")
    
  } else {
    file_name = paste(WorkingDir, "results",
                      paste0("cis_eQTLs_permutational_", cell_type, "_peer_auto"),
                      paste0("cis.permutational_", cell_type, "_FDR.significant.txt"), sep = "/")
  }
  
  significant = read.table(file_name, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = c("character", rep("NULL", 7), "character", rep("NULL", 14)))
  
  names(significant) = c("gene", "chr")
  significant = cbind(cell_type, significant)
  significant = significant[significant[, "chr"] == chr,]
  
  return(significant)
}

read_permutational = function(cell_type){
  
  if (! grepl("pos", cell_type)){
    best = overview[overview[, "cell_type"] == cell_type, "best_number_factors"]
    
    file_name = paste(WorkingDir, "results",
                      paste0("cis_eQTLs_permutational_ENSG_", cell_type, "_none_peer_", best),
                      paste0("cis.permutational_ENSG_", cell_type, "_none.txt.gz"), sep = "/")
    
  } else {
    file_name = paste(WorkingDir, "results",
                      paste0("cis_eQTLs_permutational_", cell_type, "_peer_auto"),
                      paste0("cis.permutational_", cell_type, ".txt.gz"), sep = "/")
  }
  
  permutational = read.table(file_name, header = FALSE, stringsAsFactors = FALSE, check.names = FALSE,
                             colClasses = c("character", "NULL", "integer", rep("NULL", 5), "character", "integer", rep("NULL", 11)))
  
  names(permutational) = c("gene", "transcription_start_site", "chr", "pos")
  permutational = na.omit(permutational)
  permutational = permutational[permutational[, "chr"] == chr,]
  permutational = cbind(cell_type, permutational)
  
  return(permutational)
}

compute_nominal = function(geno, expr){
  ID = colnames(geno)
  pos = as.numeric(sub(":.*", "", sub(".*?:", "", ID)))
  n = length(expr)
  my = sum(expr) / n
  m = Rfast::colmeans(geno)
  r = (Rfast::eachcol.apply(geno, expr) - n * my * m)/(n - 1)
  sx = Rfast::colVars(geno)
  beta = r / sx
  a = my - beta * m
  Se = sqrt((Rfast::colsums((t(t(geno) * beta + a) - expr)^2)/(n - 2)) / Rfast::colsums((t(t(geno) - m))^2))
  p_value = pt(abs(beta / Se), n - 2, lower.tail = FALSE) * 2
  
  nominal = data.frame(ID, pos, beta, p_value)
  
  return(nominal)
}

compute_theta = function(nominal_X, nominal_Y){
  signif_variants = nominal_X[, "p_value"] < SNP_threshold | nominal_Y[, "p_value"] < SNP_threshold
  n_variants = length(which(signif_variants))
  
  if (n_variants >= 2){
    
    if ("log_p_value" %in% names(nominal_X)){
      X_logp = nominal_X[signif_variants, "log_p_value"]
      
    } else {
      X_logp = -log10(nominal_X[signif_variants, "p_value"])
    }
    
    if ("log_p_value" %in% names(nominal_Y)){
      Y_logp = nominal_Y[signif_variants, "log_p_value"]
      
    } else {
      Y_logp = -log10(nominal_Y[signif_variants, "p_value"])
    }
    
    X_beta = nominal_X[signif_variants, "beta"]
    Y_beta = nominal_Y[signif_variants, "beta"]
    
    w_i = pmax(X_logp / max(X_logp), Y_logp / max(Y_logp))
    
    X_w = weighted.mean(X_logp, w_i)
    Y_w = weighted.mean(Y_logp, w_i)
    
    X_sigma = sqrt(sum(w_i * (X_logp - X_w) ^ 2) / sum(w_i))
    Y_sigma = sqrt(sum(w_i * (Y_logp - Y_w) ^ 2) / sum(w_i))
    
    R_w = sum(w_i * ((X_logp - X_w) / X_sigma) * ((Y_logp - Y_w) / Y_sigma)) / sum(w_i)
    R_ws = sum(w_i * ((X_logp - X_w) / X_sigma) * ((Y_logp - Y_w) / Y_sigma) * sign(X_beta * Y_beta)) / sum(w_i)
    
    theta = R_ws / (1 + exp(-30 * (R_w - 0.3)))
    
  } else {
    R_w = 0
    R_ws = 0
    theta = 0
  }
  
  if (!is.finite(theta)){
    theta = 0
  }
  
  return(c(theta, n_variants, R_w, R_ws))
}

compare_nominal = function(i, comparisons){
  print(paste(i, "over", nrow(comparisons)))
  set.seed(2023)
  
  if (analysis %in% c("EAP", "EAP_manual", "EAP_consensus")){
    gene_X = comparisons[i, "gene_X"]
    gene_Y = comparisons[i, "gene_Y"]
    left_border = comparisons[i, "left_border"]
    right_border = comparisons[i, "right_border"]
    
    if (analysis == "EAP_consensus" && paste(gene_X, cell_type_X) %in% paste(EAPs[, "gene"], EAPs[, "cell_type"])){
      module_X = EAPs[EAPs[, "gene"] == gene_X & EAPs[, "cell_type"] == cell_type_X, "module"]
      nominal_X = read.table(paste(MetaDir, paste0("Meta_analysis_for_module_", module_X, "_", chr, ".txt"), sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      nominal_X[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", nominal_X[, "ID"])))
      nominal_X = nominal_X[nominal_X[, "pos"] <= right_border & nominal_X[, "pos"] >= left_border,]
      
    } else {
      geno_X = t(genotypes_X[genotypes_X[, "pos"] <= right_border & genotypes_X[, "pos"] >= left_border, individuals_X])
      expr_X = as.numeric(expression_X[expression_X[, "gene"] == gene_X, individuals_X])
      nominal_X = compute_nominal(geno_X, expr_X)
    }
    
    if (analysis == "EAP_consensus" && paste(gene_Y, cell_type_Y) %in% paste(EAPs[, "gene"], EAPs[, "cell_type"])){
      module_Y = EAPs[EAPs[, "gene"] == gene_Y & EAPs[, "cell_type"] == cell_type_Y, "module"]
      nominal_Y = read.table(paste(MetaDir, paste0("Meta_analysis_for_module_", module_Y, "_", chr, ".txt"), sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      nominal_Y[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", nominal_Y[, "ID"])))
      nominal_Y = nominal_Y[nominal_Y[, "pos"] <= right_border & nominal_Y[, "pos"] >= left_border,]
      
    } else {
      geno_Y = t(genotypes_Y[genotypes_Y[, "pos"] <= right_border & genotypes_Y[, "pos"] >= left_border, individuals_Y])
      expr_Y = as.numeric(expression_Y[expression_Y[, "gene"] == gene_Y, individuals_Y])
      nominal_Y = compute_nominal(geno_Y, expr_Y)
    }
    
    intersection = intersect(nominal_X[, "ID"], nominal_Y[, "ID"])
    
    nominal_X = nominal_X[nominal_X[, "ID"] %in% intersection,]
    nominal_X = nominal_X[match(intersection, nominal_X[, "ID"]),]
    
    nominal_Y = nominal_Y[nominal_Y[, "ID"] %in% intersection,]
    nominal_Y = nominal_Y[match(intersection, nominal_Y[, "ID"]),]
    
    theta_n_variants = compute_theta(nominal_X, nominal_Y)
    theta = theta_n_variants[1]
    n_variants = theta_n_variants[2]
    R_w = theta_n_variants[3]
    R_ws = theta_n_variants[4]
    
    if (B == 0){
      b = B
      p_value = NA
      
    } else if (B > 0){
      early_stop = round(B * alpha) + 1
      count_X = 0
      count_Y = 0
      count = 0
      b = 0
      
      while (b < B & count < early_stop){
        perm_expr_X = as.numeric(expression_X[expression_X[, "gene"] == gene_X, sample(individuals_X)])
        perm_expr_Y = as.numeric(expression_Y[expression_Y[, "gene"] == gene_Y, sample(individuals_Y)])
        
        perm_nominal_X = compute_nominal(geno_X, perm_expr_X)
        perm_nominal_Y = compute_nominal(geno_Y, perm_expr_Y)
        
        perm_theta_X = compute_theta(perm_nominal_X, nominal_Y)[1]
        perm_theta_Y = compute_theta(nominal_X, perm_nominal_Y)[1]
        
        count_X = ifelse(abs(perm_theta_X) >= abs(theta), count_X + 1, count_X)
        count_Y = ifelse(abs(perm_theta_Y) >= abs(theta), count_Y + 1, count_Y)
        count = mean(count_X, count_Y)
        b = b + 1
      }
      
      p_value = count / b
    }
    
    result = data.frame(cell_type_X, gene_X, cell_type_Y, gene_Y, theta, chr, n_variants, b, p_value, R_w, R_ws)
    
  } else if (analysis == "DAP"){
    locus = comparisons[i, "locus"]
    gene = comparisons[i, "gene"]
    left_border = comparisons[i, "left_border"]
    right_border = comparisons[i, "right_border"]
    
    nominal_X = nominal_disease[nominal_disease[, "pos"] <= right_border &
                                  nominal_disease[, "pos"] >= left_border,]
    
    geno = t(genotypes[genotypes[, "SNPname"] %in% nominal_X[, "SNPname"], individuals])
    expr = as.numeric(expression[expression[, "gene"] == gene, individuals])
    
    nominal_Y = compute_nominal(geno, expr)
    nominal_Y[, "SNPname"] = ID_to_SNPname(nominal_Y[, "ID"])
    nominal_Y[, "beta"] = ifelse(nominal_Y[, "ID"] == nominal_Y[, "SNPname"], nominal_Y[, "beta"], nominal_Y[, "beta"] * (-1))
    
    theta_n_variants = compute_theta(nominal_X, nominal_Y)
    theta = theta_n_variants[1]
    n_variants = theta_n_variants[2]
    
    if (B == 0){
      b = B
      p_value = NA
      
    } else if (B > 0){
      early_stop = round(B * alpha) + 1
      count = 0
      b = 0
      
      while (b < B & count < early_stop){
        perm_expr = as.numeric(expression[expression[, "gene"] == gene, sample(individuals)])
        
        perm_nominal_Y = compute_nominal(geno, perm_expr)
        perm_nominal_Y[, "SNPname"] = ID_to_SNPname(perm_nominal_Y[, "ID"])
        perm_nominal_Y[, "beta"] = ifelse(perm_nominal_Y[, "ID"] == perm_nominal_Y[, "SNPname"], perm_nominal_Y[, "beta"], perm_nominal_Y[, "beta"] * (-1))
        
        perm_theta = compute_theta(nominal_X, perm_nominal_Y)[1]
        
        count = ifelse(abs(perm_theta) >= abs(theta), count + 1, count)
        b = b + 1
      }
      
      p_value = count / b
    }
    
    result = data.frame(disease, locus, cell_type, gene, theta, chr, n_variants, b, p_value)
    
  } else if (analysis == "DAP_consensus"){
    locus = comparisons[i, "locus"]
    left_border = comparisons[i, "left_border"]
    right_border = comparisons[i, "right_border"]
    
    nominal_X = nominal_disease[nominal_disease[, "pos"] <= right_border &
                                  nominal_disease[, "pos"] >= left_border,]
    
    nominal_Y = nominal_meta[nominal_meta[, "SNPname"] %in% nominal_X[, "SNPname"],]
    nominal_Y[, "beta"] = ifelse(nominal_Y[, "ID"] == nominal_Y[, "SNPname"], nominal_Y[, "beta"], nominal_Y[, "beta"] * (-1))
    
    intersection = intersect(nominal_X[, "ID"], nominal_Y[, "ID"])
    
    nominal_X = nominal_X[nominal_X[, "ID"] %in% intersection,]
    nominal_X = nominal_X[match(intersection, nominal_X[, "ID"]),]
    
    nominal_Y = nominal_Y[nominal_Y[, "ID"] %in% intersection,]
    nominal_Y = nominal_Y[match(intersection, nominal_Y[, "ID"]),]
    
    theta_n_variants = compute_theta(nominal_X, nominal_Y)
    theta = theta_n_variants[1]
    n_variants = theta_n_variants[2]
    
    result = data.frame(disease, locus, module, theta, n_variants)
  }
  
  return(result)
}

################################################################

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

args = commandArgs(trailingOnly = TRUE)
SNP_threshold = 0.05 #1.1
analysis = args[1]
dataset = args[2]

if (dataset == "blood"){
  WorkingDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA", sep = "/")
  
} else if (dataset == "biopsies"){
  WorkingDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap", sep = "/")
  
} else if (dataset == "both"){
  WorkingDir = c(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/nextflow_eQTL_pipeline_12e_PCA", sep = "/"),
                 paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/CROHN/Viacheslav/Cedar2/Transcriptome/scData_analysis/cis_eap_dap", sep = "/"))
}

ThetaDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta", dataset, sep = "/")
DiseaseDir = paste(prefix, "Analysis/Theta_correlation/Pipeline_test/docs", sep = "/")
directories = list.files(path = paste(WorkingDir, "results", sep = "/"), pattern = "cis_eQTLs_permutational_", include.dirs = TRUE)
directories = gsub("cis_eQTLs_permutational_|ENSG_|none_", "", directories)

overview = data.frame(cell_type = sub("_peer_.*", "", directories),
                      best_number_factors = sub(".*_peer_", "", directories))

if (analysis %in% c("EAP", "EAP_manual", "EAP_consensus")){
  cell_type_X = args[3]
  cell_type_Y = args[4]
  chr = args[5]
  B = as.numeric(args[6])
  alpha = as.numeric(args[7])
  
  EAP_sym_border = 10^6
  cell_type = unique(c(cell_type_X, cell_type_Y))
  
  permutational = do.call("rbind", lapply(cell_type, read_permutational))
  title = paste0("theta_", paste(cell_type, collapse = "_"), "_", chr, ".txt")
  OutputDir = paste(ThetaDir, analysis, paste0(B, "_permutations", ifelse(B != 0, paste0("_", alpha), "")), cell_type_X, sep = "/")
  
} else if (analysis == "DAP"){
  disease = args[3]
  cell_type = args[4]
  chr = args[5]
  B = as.numeric(args[6])
  alpha = as.numeric(args[7])
  DAP_borders = args[8]
  DAP_sym_border = as.numeric(args[9])
  
  permutational = do.call("rbind", lapply(cell_type, read_permutational))
  loci = readRDS(paste(DiseaseDir, paste0("loci_", disease, "_", DAP_borders, ".rds"), sep = "/"))
  loci[, "left_border"] = loci[, "left_border"] - DAP_sym_border
  loci[, "right_border"] = loci[, "right_border"] + DAP_sym_border
  title = paste0("theta_", disease, "_", cell_type, "_", chr, ".txt")
  OutputDir = paste(ThetaDir, analysis, paste(B, "permutations", alpha, DAP_borders, format(DAP_sym_border, scientific = FALSE), sep = "_"), cell_type, sep = "/")
  
} else if (analysis == "DAP_consensus"){
  disease = args[3]
  module = args[4]
  chr = args[5]
  DAP_borders = args[6]
  DAP_sym_border = as.numeric(args[7])
  iteration = args[8]
  
  loci = readRDS(paste(DiseaseDir, paste0("loci_", disease, "_", DAP_borders, ".rds"), sep = "/"))
  loci[, "left_border"] = loci[, "left_border"] - DAP_sym_border
  loci[, "right_border"] = loci[, "right_border"] + DAP_sym_border
  
  MetaDir = paste(ThetaDir, paste0("meta_analysis_regulatory_", iteration, "_extended"), sep = "/")
  nominal_meta = read.table(paste(MetaDir, paste0("Meta_analysis_for_module_", module, "_", chr, ".txt"), sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  nominal_meta[, "SNPname"] = ID_to_SNPname(nominal_meta[, "ID"])
  nominal_meta[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", nominal_meta[, "ID"])))
  nominal_meta = nominal_meta[order(nominal_meta[, "pos"]),]
  
  title = paste0("theta_", disease, "_", module, "_", chr, ".txt")
  OutputDir = paste(ThetaDir, analysis, paste(DAP_borders, format(DAP_sym_border, scientific = FALSE), sep = "_"), sep = "/")
}

if (!dir.exists(OutputDir)){
  dir.create(OutputDir, recursive = TRUE)
}

if (analysis == "EAP_consensus"){
  MetaDir = paste(ThetaDir, paste0("meta_analysis_gene_specific_extended"), sep = "/")
  EAPs = readRDS(paste(ThetaDir, paste0(dataset, "_EAPs_representative_gene_specific.rds"), sep = "/"))
  EAPs = EAPs[EAPs[, "cell_type"] %in% cell_type & EAPs[, "chromosome_name"] == sub("chr", "", chr) & EAPs[, "nr_EAPs"] > 1, c("cell_type", "gene", "module")]
}

################################################################

comparisons = NULL

if (analysis %in% c("EAP", "EAP_consensus")){
  significant = do.call("rbind", lapply(cell_type, read_significant))
  
  if (nrow(significant) != 0){
    permutational[, "left_border"] = permutational[, "transcription_start_site"] - EAP_sym_border
    permutational[, "right_border"] = permutational[, "transcription_start_site"] + EAP_sym_border
    
    for (gene in significant[significant[, "cell_type"] == cell_type_X, "gene"]){
      left_border = permutational[permutational[, "cell_type"] == cell_type_X & permutational[, "gene"] == gene, "left_border"]
      right_border = permutational[permutational[, "cell_type"] == cell_type_X & permutational[, "gene"] == gene, "right_border"]
      
      close_genes = permutational[permutational[, "cell_type"] == cell_type_Y &
                                    permutational[, "left_border"] <= right_border &
                                    permutational[, "right_border"] >= left_border, c("gene", "left_border", "right_border")]
      
      if (nrow(close_genes) > 0){
        comparisons = rbind(comparisons, data.frame(gene_X = gene, gene_Y = close_genes[, "gene"],
                                                    left_border = pmin(left_border, close_genes[, "left_border"]),
                                                    right_border = pmax(right_border, close_genes[, "right_border"])))
      }
    }
    
    if (cell_type_X == cell_type_Y){
      alpha_order = paste(do.call(pmin, comparisons[, c("gene_X", "gene_Y")]),
                          do.call(pmax, comparisons[, c("gene_X", "gene_Y")]), sep = "_")
      
      comparisons = comparisons[!duplicated(alpha_order),]
      comparisons = comparisons[comparisons[, "gene_X"] != comparisons[, "gene_Y"],]
      
    } else {
      
      for (gene in significant[significant[, "cell_type"] == cell_type_Y, "gene"]){
        left_border = permutational[permutational[, "cell_type"] == cell_type_Y & permutational[, "gene"] == gene, "left_border"]
        right_border = permutational[permutational[, "cell_type"] == cell_type_Y & permutational[, "gene"] == gene, "right_border"]
        
        close_genes = permutational[permutational[, "cell_type"] == cell_type_X &
                                      permutational[, "left_border"] <= right_border &
                                      permutational[, "right_border"] >= left_border, c("gene", "left_border", "right_border")]
        
        if (nrow(close_genes) > 0){
          comparisons = rbind(comparisons, data.frame(gene_X = close_genes[, "gene"], gene_Y = gene,
                                                      left_border = pmin(left_border, close_genes[, "left_border"]),
                                                      right_border = pmax(right_border, close_genes[, "right_border"])))
        }
      }
      
      comparisons = unique(comparisons)
    }
  }
  
} else if (analysis == "EAP_manual"){
  EAPs = readRDS(paste(ThetaDir, paste0(dataset, "_EAPs_representative_regulatory_first.rds"), sep = "/"))
  EAPs[, "left_border"] = EAPs[, "transcription_start_site"] - EAP_sym_border
  EAPs[, "right_border"] = EAPs[, "transcription_start_site"] + EAP_sym_border
  EAPs = EAPs[EAPs[, "chromosome_name"] == sub("chr", "", chr),]
  EAPs = EAPs[EAPs[, "signif"] == FALSE,]
  EAPs = EAPs[EAPs[, "cell_type"] %in% cell_type,]
  
  if (nrow(EAPs) > 0){
    
    for (gene in EAPs[EAPs[, "cell_type"] == cell_type_X, "gene"]){
      module = EAPs[EAPs[, "cell_type"] == cell_type_X & EAPs[, "gene"] == gene, "module"]
      left_border = EAPs[EAPs[, "cell_type"] == cell_type_X & EAPs[, "gene"] == gene, "left_border"]
      right_border = EAPs[EAPs[, "cell_type"] == cell_type_X & EAPs[, "gene"] == gene, "right_border"]
      
      close_genes = EAPs[EAPs[, "cell_type"] == cell_type_Y &
                           EAPs[, "module"] == module &
                           EAPs[, "left_border"] <= right_border &
                           EAPs[, "right_border"] >= left_border, c("gene", "left_border", "right_border")]
      
      if (nrow(close_genes) > 0){
        comparisons = rbind(comparisons, data.frame(gene_X = gene, gene_Y = close_genes[, "gene"],
                                                    left_border = pmin(left_border, close_genes[, "left_border"]),
                                                    right_border = pmax(right_border, close_genes[, "right_border"])))
      }
    }
    
    if (cell_type_X == cell_type_Y){
      alpha_order = paste(do.call(pmin, comparisons[, c("gene_X", "gene_Y")]),
                          do.call(pmax, comparisons[, c("gene_X", "gene_Y")]), sep = "_")
      
      comparisons = comparisons[!duplicated(alpha_order),]
      comparisons = comparisons[comparisons[, "gene_X"] != comparisons[, "gene_Y"],]
    }
  }
  
} else if (analysis == "DAP"){
  
  for (locus in loci[loci[, "chr"] == chr, "rsID"]){
    left_border = loci[loci[, "rsID"] == locus, "left_border"]
    right_border = loci[loci[, "rsID"] == locus, "right_border"]
    
    close_genes = permutational[permutational[, "pos"] <= right_border &
                                  permutational[, "pos"] >= left_border, "gene"]
    
    if (length(close_genes) > 0){
      comparisons = rbind(comparisons, data.frame(locus, gene = close_genes, left_border, right_border))
    }
  }
  
} else if (analysis == "DAP_consensus"){
  top_SNP_pos = nominal_meta[which.max(nominal_meta[, "abs_z_score"]), "pos"]
  
  loci = loci[loci[, "chr"] == chr & loci[, "left_border"] <= top_SNP_pos & loci[, "right_border"] >= top_SNP_pos,]
  
  if (nrow(loci) > 0){
    comparisons = data.frame(locus = loci[, "rsID"], module, left_border = loci[, "left_border"], right_border = loci[, "right_border"])
  }
}

if (!is.null(comparisons)){
  
  if (nrow(comparisons) > 0){
    
    if (! file.exists(paste(OutputDir, title, sep = "/")) || file.size(paste(OutputDir, title, sep = "/")) == 0 || nrow(comparisons) != (length(readLines(paste(OutputDir, title, sep = "/"))) - 1)){
      
      if (analysis %in% c("EAP", "EAP_manual", "EAP_consensus")){
        expression_X = read_expression(cell_type_X)
        expression_Y = read_expression(cell_type_Y)
        individuals_X = setdiff(names(expression_X), c("cell_type", "#chr", "gene"))
        individuals_Y = setdiff(names(expression_Y), c("cell_type", "#chr", "gene"))
        genotypes_X = read_genotypes(cell_type_X, individuals_X)
        genotypes_Y = read_genotypes(cell_type_Y, individuals_Y)
        genotypes_X = genotypes_X[row.names(genotypes_X) %in% row.names(genotypes_Y),]
        genotypes_Y = genotypes_Y[row.names(genotypes_Y) %in% row.names(genotypes_X),]
        variants_X = row.names(genotypes_X[genotypes_X[, "MAF"] >= 0.05 & genotypes_X[, "MAF"] <= 0.95,])
        variants_Y = row.names(genotypes_Y[genotypes_Y[, "MAF"] >= 0.05 & genotypes_Y[, "MAF"] <= 0.95,])
        genotypes_X = genotypes_X[row.names(genotypes_X) %in% union(variants_X, variants_Y),]
        genotypes_Y = genotypes_Y[row.names(genotypes_Y) %in% union(variants_Y, variants_X),]
        
      } else if (analysis == "DAP"){
        expression = read_expression(cell_type)
        individuals = setdiff(names(expression), c("cell_type", "#chr", "gene"))
        genotypes = read_genotypes(cell_type, individuals)
        genotypes = genotypes[genotypes[, "MAF"] >= 0.05 & genotypes[, "MAF"] <= 0.95,]
        genotypes[, "SNPname"] = ID_to_SNPname(row.names(genotypes))
        nominal_disease = readRDS(paste(DiseaseDir, paste0("nominal_raw_", disease, "_subset.rds"), sep = "/"))
        nominal_disease = nominal_disease[nominal_disease[, "SNPname"] %in% genotypes[, "SNPname"],]
        
      } else if (analysis == "DAP_consensus"){
        nominal_disease = readRDS(paste(DiseaseDir, paste0("nominal_raw_", disease, "_subset.rds"), sep = "/"))
        nominal_disease = nominal_disease[nominal_disease[, "SNPname"] %in% nominal_meta[, "SNPname"],]
      }
      
      results = do.call("rbind", lapply(1:nrow(comparisons), compare_nominal, comparisons))
      write.table(results, paste(OutputDir, title, sep = "/"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
  }
}
