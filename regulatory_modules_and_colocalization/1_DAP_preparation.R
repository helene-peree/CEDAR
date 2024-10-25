library(rtracklayer)
library(data.table)

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

ID_to_SNPname = function(ID){
  coordinates = data.frame(do.call("rbind", strsplit(ID, ":")))
  names(coordinates) = c("chr", "pos", "ref", "alt")
  
  SNPname = paste(coordinates[, "chr"], coordinates[, "pos"],
                  do.call(pmin, coordinates[, c("ref", "alt")]),
                  do.call(pmax, coordinates[, c("ref", "alt")]), sep = ":")
  
  return(SNPname)
}

#############

if (grepl("Y:", getwd())){
  prefix = "Y:" # RStudio
  
} else {
  prefix = "~" # GIGA cluster
}

summary_dir = paste(prefix, "_SHARE_/Research/GEN/UAG/IBD/PLATFORMS/GEN/Slava_pipeline/gwas_summary", sep = "/")
output_file_dir = paste(prefix, "Analysis/Theta_correlation/Pipeline_test/docs", sep = "/")

if (!dir.exists(output_file_dir)) {
  dir.create(output_file_dir, recursive = TRUE)
}

CD_input_file = paste(summary_dir, "cd_build37_40266_20161107.txt.gz", sep = "/")
UC_input_file = paste(summary_dir, "uc_build37_45975_20161107.txt.gz", sep = "/")
IBD_input_file = paste(summary_dir, "ibd_build37_59957_20161107.txt.gz", sep = "/")

CD_output_file = paste(output_file_dir, "nominal_raw_CD_subset.rds", sep = "/")
UC_output_file = paste(output_file_dir, "nominal_raw_UC_subset.rds", sep = "/")
IBD_output_file = paste(output_file_dir, "nominal_raw_IBD_subset.rds", sep = "/")

chain_file = paste(summary_dir, "hg19ToHg38.over.chain", sep = "/")
chain_info = rtracklayer::import.chain(chain_file)

for (disease in c("CD", "UC", "IBD")){
  
  if (disease == "CD"){
    nominal = data.table::fread(CD_input_file)
    
  } else if (disease == "UC"){
    nominal = data.table::fread(UC_input_file)
    
  } else if (disease == "IBD"){
    nominal = data.table::fread(IBD_input_file)
  }
  
  nominal = as.data.frame(nominal, stringsAsFactors = FALSE)
  nominal[, "chr"] = sub(":.*", "", nominal[, "MarkerName"])
  nominal[, "pos"] = as.numeric(sub(".*:", "", sub("_.*", "", nominal[, "MarkerName"])))
  nominal = lift_over(nominal, "pos")
  nominal[, "ID"] = paste(nominal[, "chr"], nominal[, "pos"], toupper(nominal[, "Allele1"]), toupper(nominal[, "Allele2"]), sep = ":")
  names(nominal)[names(nominal) == "Effect"] = "beta"
  names(nominal)[names(nominal) == "P.value"] = "p_value"
  nominal[, "SNPname"] = ID_to_SNPname(nominal[, "ID"])
  nominal[, "beta"] = ifelse(nominal[, "ID"] == nominal[, "SNPname"], nominal[, "beta"], nominal[, "beta"] * (-1))
  nominal = nominal[, c("chr", "pos", "ID", "SNPname", "beta", "p_value")]
  chr_numbers = as.numeric(sub("chr", "", sub("X", "23", nominal[, "chr"])))
  nominal = nominal[order(chr_numbers, nominal[, "pos"]),]
  
  saveRDS(nominal, paste(output_file_dir, paste0("nominal_raw_", disease, "_subset.rds"), sep = "/"))
}
