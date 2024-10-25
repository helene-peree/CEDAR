args = commandArgs(trailingOnly = TRUE)
dataset = args[1] # "blood" or "biopsies" or "both"
analysis = args[2] # "gene_specific" or "regulatory"
iteration = args[3] # "first" or "second"
module = args[4]
extended = FALSE # TRUE or FALSE

################################################################

suppressMessages(library(Rfast))

################################################################

source("function_definitions.R")

compute_nominal = function(geno, expr){
  ID = colnames(geno)
  n = length(expr)
  my = sum(expr) / n
  mx = Rfast::colmeans(geno)
  r = (Rfast::eachcol.apply(geno, expr) - n * my * mx)/(n - 1)
  sx = Rfast::colVars(geno)
  beta = r / sx
  alpha = my - beta * mx
  se = sqrt((Rfast::colsums((t(t(geno) * beta + alpha) - expr)^2)/(n - 2)) / Rfast::colsums((t(t(geno) - mx))^2))
  p_value = 2 * pt(abs(beta / se), n - 2, lower.tail = FALSE)
  log_p_value = -log10(p_value)
  
  nominal = data.frame(ID, n, beta, se, p_value, log_p_value)
  
  return(nominal)
}

################################################################

MetaDir = paste(ThetaDir, paste0("meta_analysis_", analysis, iter, ifelse(extended, paste0("_extended"), "")), sep = "/")
OutputDir = paste(prefix, "Analysis/eQTL/scripts/Website", paste0(dataset, "_meta_analysis_", analysis, iter, ifelse(extended, paste0("_extended"), "")), sep = "/")

if (!dir.exists(OutputDir)){
  dir.create(OutputDir, recursive = TRUE)
}

if (!dir.exists(MetaDir)){
  dir.create(MetaDir, recursive = TRUE)
}

EAPs = readRDS(paste(ThetaDir, paste0(dataset, "_EAPs_representative_", analysis, iter, ".rds"), sep = "/"))
EAPs = EAPs[EAPs[, "module"] == module,]
EAPs[, "dataset"] = ifelse(! grepl("_pos_", EAPs[, "cell_type"]), "blood", "biopsies")

################################################################

if (nrow(EAPs) > 1){
  
  if (extended){
    distance = 3 * 10^6
    
  } else {
    distance = 10^6
  }
  
  left_border = min(EAPs[, "transcription_start_site"]) - distance
  right_border = max(EAPs[, "transcription_start_site"]) + distance
  chr = paste0("chr", unique(EAPs[, "chromosome_name"]))
  
  nominal = NULL
  
  for (new_dataset in unique(EAPs[, "dataset"])){
    genotypes = read_genotypes(new_dataset)
    genotypes[, "pos"] = as.numeric(sub(":.*", "", sub(".*?:", "", row.names(genotypes))))
    genotypes = genotypes[genotypes[, "pos"] >= left_border & genotypes[, "pos"] <= right_border,]
    
    for (cell_type in unique(EAPs[EAPs[, "dataset"] == new_dataset, "cell_type"])){
      EAPs_subset = EAPs[EAPs[, "dataset"] == new_dataset & EAPs[, "cell_type"] == cell_type,]
      expression = read_expression(new_dataset, cell_type, "transformed")
      individuals = setdiff(names(expression), c("cell_type", "#chr", "gene"))
      
      geno = genotypes
      geno[, "MAF"] = rowSums(round(geno[, individuals])) / (2 * length(individuals))
      variants = row.names(geno[geno[, "MAF"] >= 0.05 & geno[, "MAF"] <= 0.95,])
      geno = t(geno[, individuals])
      
      for (gene in unique(EAPs_subset[, "gene"])){
        known = EAPs_subset[EAPs_subset[, "gene"] == gene, "known"]
        expr = as.numeric(expression[expression[, "gene"] == gene, individuals])
        nominal_temp = cbind(compute_nominal(geno, expr), cell_type, gene)
        
        if (! is.na(known) & known == "representative"){
          representative = nominal_temp
        }
        
        nominal_temp = nominal_temp[nominal_temp[, "ID"] %in% variants,]
        nominal = rbind(nominal, nominal_temp)
      }
    }
  }
  
  if (dataset == "both"){
    nominal = nominal[nominal[, "ID"] %in% representative[, "ID"],]
  }
  
  nominal[, "z"] = qnorm(nominal[, "p_value"] / 2, lower.tail = FALSE)
  
  nominal_meta = NULL
  
  for (ID in sort(unique(nominal[, "ID"]))){
    beta = sign(representative[representative[, "ID"] == ID, "beta"])
    z_score = sum(nominal[nominal[, "ID"] == ID, "z"]^2)
    abs_z_score = abs(z_score)
    p_value = pchisq(abs_z_score, df = nrow(EAPs), ncp = 0, lower.tail = FALSE)
    log_p_value = -log10(p_value)
    
    nominal_meta = rbind(nominal_meta, data.frame(ID, beta, z_score, abs_z_score, p_value, log_p_value))
  }
  
  if (any(!is.finite(nominal_meta[, "log_p_value"]))){
    print("predict p_values")
    nominal_meta = nominal_meta[order(abs(nominal_meta[, "abs_z_score"])),]
    lm_results = loess(nominal_meta[is.finite(nominal_meta[, "log_p_value"]), "log_p_value"] ~ nominal_meta[is.finite(nominal_meta[, "log_p_value"]), "abs_z_score"], control = loess.control(surface = "direct"))
    x = nominal_meta[is.finite(nominal_meta[, "log_p_value"]), "abs_z_score"]
    y = lm_results$fitted
    nominal_meta[!is.finite(nominal_meta[, "log_p_value"]), "log_p_value"] = predict(lm_results, nominal_meta[!is.finite(nominal_meta[, "log_p_value"]), "abs_z_score"], se = TRUE)$fit
  }
  
  filename = paste0("Meta_analysis_for_module_", module, "_", chr, ".txt")
  write.table(nominal_meta, paste(MetaDir, filename, sep = "/"), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
