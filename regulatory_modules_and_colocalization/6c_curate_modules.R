args = commandArgs(trailingOnly = TRUE)
dataset = args[1] # "blood" or "biopsies" or "both"
analysis = args[2] # "gene_specific" or "regulatory"
iteration = args[3] # "first" or "second"
chr = args[4] # from "chr1" to "chrX"

################################################################

suppressMessages(library(Rfast))

################################################################

source("function_definitions.R")

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
    X_logp = -log10(nominal_X[signif_variants, "p_value"])
    Y_logp = nominal_Y[signif_variants, "log_p_value"]
    
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

################################################################

SNP_threshold = 0.05

OutputDir = paste(ThetaDir, paste0("check_", analysis, iter), sep = "/")

if (!dir.exists(OutputDir)){
  dir.create(OutputDir, recursive = TRUE)
}

EAPs = readRDS(paste(ThetaDir, paste0(dataset, "_EAPs_representative_", analysis, iter, ".rds"), sep = "/"))
EAPs = EAPs[EAPs[, "nr_EAPs"] > 1,]
EAPs[, "dataset"] = ifelse(! grepl("_pos_", EAPs[, "cell_type"]), "blood", "biopsies")

for (new_dataset in unique(EAPs[, "dataset"])){
  print(new_dataset)
  genotypes = read_genotypes(new_dataset)
  genotypes[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", row.names(genotypes))))
  
  for (cell_type in unique(EAPs[EAPs[, "dataset"] == new_dataset & EAPs[, "chromosome_name"] == sub("chr", "", chr), "cell_type"])){
    expression = read_expression(new_dataset, cell_type, "transformed")
    individuals = setdiff(names(expression), c("cell_type", "#chr", "gene"))
    results = NULL
    
    for (module in unique(EAPs[EAPs[, "dataset"] == new_dataset & EAPs[, "chromosome_name"] == sub("chr", "", chr) & EAPs[, "cell_type"] == cell_type, "module"])){
      print(module)
      filename = paste0("Meta_analysis_for_module_", module, "_", chr, ".txt")
      nominal_meta = read.table(paste(MetaDir, filename, sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
      nominal_meta[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", nominal_meta[, "ID"])))
      nominal_meta = nominal_meta[order(nominal_meta[, "pos"]),]
      
      left_border = min(nominal_meta[, "pos"])
      right_border = max(nominal_meta[, "pos"])
      
      geno = genotypes
      geno[, "MAF"] = rowSums(round(geno[, individuals])) / (2 * length(individuals))
      geno = geno[geno[, "MAF"] >= 0.05 & geno[, "MAF"] <= 0.95,]
      geno = t(geno[geno[, "pos"] <= right_border & geno[, "pos"] >= left_border, individuals])
      
      for (gene in EAPs[EAPs[, "dataset"] == new_dataset & EAPs[, "module"] == module & EAPs[, "cell_type"] == cell_type, "gene"]){
        expr = as.numeric(expression[expression[, "gene"] == gene, individuals])
        nominal = compute_nominal(geno, expr)
        common_variants = intersect(nominal[, "ID"], nominal_meta[, "ID"])
        
        theta_n_variants = compute_theta(nominal[nominal[, "ID"] %in% common_variants,], nominal_meta[nominal_meta[, "ID"] %in% common_variants,])
        theta = theta_n_variants[1]
        n_variants = theta_n_variants[2]
        R_w = theta_n_variants[3]
        R_ws = theta_n_variants[4]
        
        results = rbind(results, data.frame(gene, cell_type, chr, theta, n_variants, R_w, R_ws))
      }
    }
    
    write.table(results, paste(OutputDir, paste0("theta_", chr, "_", cell_type, ".txt"), sep = "/"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  }
}
