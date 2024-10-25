suppressMessages(library(Rfast))
library(rtracklayer)

################################################################

ID_to_SNPname = function(ID){
  coordinates = data.frame(do.call("rbind", strsplit(ID, ":")))
  names(coordinates) = c("chr", "pos", "ref", "alt")
  
  SNPname = paste(coordinates[, "chr"], coordinates[, "pos"],
                  do.call(pmin, coordinates[, c("ref", "alt")]),
                  do.call(pmax, coordinates[, c("ref", "alt")]), sep = ":")
  
  return(SNPname)
}

lift_over = function(dataset, column){
  fake_column = paste0("fake_", column)
  excluded_columns = c(fake_column, "end", "width", "strand")
  dataset[, fake_column] = dataset[, column] + 1
  dataset_GR = GenomicRanges::makeGRangesFromDataFrame(dataset,
                                                       keep.extra.columns = T,
                                                       ignore.strand = F,
                                                       seqnames.field = "chr",
                                                       start.field = column,
                                                       end.field = fake_column)
  seqlevelsStyle(dataset_GR) = "UCSC"
  dataset = rtracklayer::liftOver(dataset_GR, chain_info)
  dataset = as.data.frame(unlist(dataset), stringsAsFactors = FALSE)
  dataset = dataset[, setdiff(names(dataset), excluded_columns)]
  names(dataset)[names(dataset) == "seqnames"] = "chr"
  names(dataset)[names(dataset) == "start"] = column
  
  return(dataset)
}

read_Franke = function(){
  
  if (! file.exists(paste(FrankeDir, paste0("2019-12-11-cis-eQTLsFDR-ProbeLevel-GRCh38_", chr, ".rds"), sep = "/"))){
    file_name = "2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
    
    dataset = read.table(paste(FrankeDir, file_name, sep = "/"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    names(dataset)[names(dataset) == "Pvalue"] = "p_value"
    names(dataset)[names(dataset) == "SNPChr"] = "chr"
    names(dataset)[names(dataset) == "SNPPos"] = "pos"
    names(dataset)[names(dataset) == "AssessedAllele"] = "ref"
    names(dataset)[names(dataset) == "OtherAllele"] = "alt"
    names(dataset)[names(dataset) == "Zscore"] = "z_score"
    names(dataset)[names(dataset) == "Gene"] = "gene"
    names(dataset)[names(dataset) == "GeneSymbol"] = "gene_name"
    names(dataset)[names(dataset) == "GenePos"] = "transcription_start_site"
    dataset = dataset[, c("p_value", "chr", "pos", "ref", "alt", "z_score", "gene", "gene_name", "transcription_start_site", "NrCohorts", "FDR")]
    
    dataset = lift_over(dataset, "pos")
    dataset = dataset[dataset[, "chr"] == chr,]
    dataset[, "abs_z_score"] = abs(dataset[, "z_score"])
    dataset[, "p_value"] = 2 * pnorm(dataset[, "abs_z_score"], lower.tail = FALSE)
    dataset[, "log_p_value"] = -log10(dataset[, "p_value"])
    
    if (any(!is.finite(dataset[, "log_p_value"]))){
      print("predict p_values")
      dataset = dataset[order(dataset[, "abs_z_score"]),]
      dataset_subset = dataset[is.finite(dataset[, "log_p_value"]),]
      nrow(dataset_subset)
      print(seq(1, nrow(dataset_subset), by = nrow(dataset_subset)/1000))
      dataset_subset = dataset_subset[seq(1, nrow(dataset_subset), by = nrow(dataset_subset)/1000),]
      nrow(dataset_subset)
      
      lm_results = loess(dataset_subset[, "log_p_value"] ~ dataset_subset[, "abs_z_score"], control = loess.control(surface = "direct"))
      x = dataset_subset[, "abs_z_score"]
      y = lm_results$fitted
      dataset[!is.finite(dataset[, "log_p_value"]), "log_p_value"] = predict(lm_results, dataset[!is.finite(dataset[, "log_p_value"]), "abs_z_score"], se = TRUE)$fit
    }
    
    saveRDS(dataset, paste(FrankeDir, paste0("2019-12-11-cis-eQTLsFDR-ProbeLevel-GRCh38_", chr, ".rds"), sep = "/"))
    
  } else {
    dataset = readRDS(paste(FrankeDir, paste0("2019-12-11-cis-eQTLsFDR-ProbeLevel-GRCh38_", chr, ".rds"), sep = "/"))
  }
  
  return(dataset)
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
    Y_beta = nominal_Y[signif_variants, "z_score"]
    
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
  
  if (analysis == "DAP"){
    locus = comparisons[i, "locus"]
    gene = comparisons[i, "gene"]
    left_border = comparisons[i, "left_border"]
    right_border = comparisons[i, "right_border"]
    
    nominal_Y = Franke[Franke[, "gene"] == gene & Franke[, "pos"] <= right_border & Franke[, "pos"] >= left_border,]
    nominal_X = nominal_disease[nominal_disease[, "SNPname"] %in% nominal_Y[, "SNPname"],]
    
    theta_n_variants = compute_theta(nominal_X, nominal_Y)
    theta = theta_n_variants[1]
    n_variants = theta_n_variants[2]
    
    b = 0
    p_value = NA
    
    result = data.frame(disease, locus, cell_type, gene, theta, chr, n_variants, b, p_value)
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

FrankeDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/eQTL/Franke", sep = "/")
WorkingDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta/Franke", sep = "/")
DiseaseDir = paste(prefix, "Analysis/Theta_correlation/Pipeline_test/docs", sep = "/")

chain_info = rtracklayer::import.chain(paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/PLATFORMS/GEN/Slava_pipeline/gwas_summary/hg19ToHg38.over.chain", sep = "/"))

if (analysis == "DAP"){
  disease = args[2]
  cell_type = "Franke"
  chr = args[3]
  DAP_borders = args[4]
  DAP_sym_border = as.numeric(args[5])
  
  Franke = read_Franke()
  Franke = Franke[order(Franke[, "pos"]),]
  Franke[, "ID"] = paste(Franke[, "chr"], Franke[, "pos"], Franke[, "ref"], Franke[, "alt"], sep = ":")
  Franke[, "SNPname"] = ID_to_SNPname(Franke[, "ID"])
  Franke[, "z_score"] = ifelse(Franke[, "ID"] == Franke[, "SNPname"], Franke[, "z_score"], Franke[, "z_score"] * (-1))
  
  loci = readRDS(paste(DiseaseDir, paste0("loci_", disease, "_", DAP_borders, ".rds"), sep = "/"))
  title = paste0("theta_", disease, "_", cell_type, "_", chr, ".txt")
  OutputDir = paste(WorkingDir, analysis, paste(DAP_borders, format(DAP_sym_border, scientific = FALSE), sep = "_"), sep = "/")
}

if (!dir.exists(OutputDir)){
  dir.create(OutputDir, recursive = TRUE)
}

permutational = data.frame(cell_type, unique(Franke[, c("gene", "chr")]))
permutational[, "pos"] = sapply(permutational[, "gene"], function(x){
  x = Franke[Franke[, "gene"] == x, c("abs_z_score", "pos")]
  x[which.max(x[, "abs_z_score"]), "pos"]})

significant = data.frame(cell_type, unique(Franke[Franke[, "FDR"] < 0.05, c("gene", "chr")]))

################################################################

if (! file.exists(paste(OutputDir, title, sep = "/"))){
  comparisons = NULL
  
  if (analysis == "DAP"){
    
    for (locus in loci[loci[, "chr"] == chr, "rsID"]){
      left_border = loci[loci[, "rsID"] == locus, "left_border"] - DAP_sym_border
      right_border = loci[loci[, "rsID"] == locus, "right_border"] + DAP_sym_border
      
      close_genes = permutational[permutational[, "pos"] <= right_border &
                                    permutational[, "pos"] >= left_border, "gene"]
      
      if (length(close_genes) > 0){
        comparisons = rbind(comparisons, data.frame(locus, gene = close_genes, left_border, right_border))
      }
    }
  }
  
  if (!is.null(comparisons)){
    
    if (nrow(comparisons) > 0){
      
      if (analysis == "DAP"){
        nominal_disease = readRDS(paste(DiseaseDir, paste0("nominal_raw_", disease, "_subset.rds"), sep = "/"))
        nominal_disease = nominal_disease[nominal_disease[, "SNPname"] %in% Franke[, "SNPname"],]
        Franke = Franke[Franke[, "SNPname"] %in% nominal_disease[, "SNPname"],]
      }
      
      results = do.call("rbind", lapply(1:nrow(comparisons), compare_nominal, comparisons))
      write.table(results, paste(OutputDir, title, sep = "/"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
    }
  }
}
