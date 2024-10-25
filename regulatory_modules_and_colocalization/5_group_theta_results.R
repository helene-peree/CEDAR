args = commandArgs(trailingOnly = TRUE)
dataset = args[1] # "blood" or "biopsies" or "both"
analysis = args[2] # "EAP" or "EAP_manual" or "EAP_consensus" or "DAP" or "DAP_consensus"

##################################

concatenate_chrs = function(comparison){
  comparison_files = grep(comparison, files, value = TRUE)
  
  if (length(comparison_files) > 0){
    
    if (all(file.size(comparison_files) > 0)){
      print(basename(comparison))
      comparison_data = do.call("rbind", lapply(comparison_files, read.table, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE))
      comparison_data = comparison_data[comparison_data[, "n_variants"] >= 2,]
      write.table(comparison_data, sub("_chr", ".tsv", comparison), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      file.remove(comparison_files)
      
    } else {
      print(paste(basename(comparison_files[file.size(comparison_files) == 0]), "problem file size"))
      print(error)
    }
  }
}

##################################

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

WorkingDir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/SYSCID/Helene/theta", dataset, sep = "/")

if (analysis == "EAP"){
  InputDir = paste(WorkingDir, analysis, "0_permutations", sep = "/")
  pattern = ".tsv"
  
  if (dataset == "blood"){
    filename = paste0(dataset, "_thetas_", analysis, ".rds")
    
  } else if (dataset %in% c("biopsies", "both")){
    filename = paste0(dataset, "_thetas_", analysis, "_above_threshold.rds")
  }
  
} else if (analysis %in% c("EAP_manual", "EAP_consensus")){
  InputDir = paste(WorkingDir, analysis, "0_permutations", sep = "/")
  pattern = ".tsv"
  filename = paste0(dataset, "_thetas_", analysis, ".rds")
  
} else if (analysis == "DAP"){
  InputDir = paste(WorkingDir, analysis, "10000_permutations_0.1_manual_0", sep = "/")
  pattern = ".tsv"
  filename = paste0(dataset, "_thetas_", analysis, ".rds")
  
} else if (analysis == "DAP_consensus"){
  InputDir = paste(WorkingDir, analysis, "manual_0", sep = "/")
  pattern = ".txt"
  filename = paste0(dataset, "_thetas_", analysis, ".rds")
}

##################################

if (analysis != "DAP_consensus"){
  directories = list.dirs(path = InputDir, full.names = TRUE, recursive = FALSE)
  
  for (directory in directories){
    files = list.files(path = directory, pattern = ".txt", full.names = TRUE, recursive = FALSE)
    
    if (length(files) > 0){
      lapply(unique(sub("chr.*", "chr", files)), concatenate_chrs)
    }
  }
}

##################################

files = list.files(path = InputDir, pattern = pattern, full.names = TRUE, recursive = TRUE)

coloc = do.call("rbind", lapply(files, function(file_name){
  x = read.table(file_name, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  if (grepl("above_threshold", filename)){
    x = x[abs(x[, "theta"]) >= 0.6,]
  }
}))

saveRDS(coloc, paste(WorkingDir, filename, sep = "/"))
