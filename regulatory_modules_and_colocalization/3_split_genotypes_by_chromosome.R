read_genotypes = function(RawDir, RawFile){
  file_name = paste(RawDir, paste0(RawFile, ".vcf.gz"), sep = "/")
  
  vcf = file(file_name, "r")
  skip = 0
  line = readLines(vcf, 1)
  
  while(!grepl("#CHROM", line)){
    skip = skip + 1
    line = readLines(vcf, 1)
  }
  
  close(vcf)
  
  genotypes = read.table(file_name, comment.char = "", skip = skip, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(genotypes) = genotypes[, "ID"]
  genotypes = as.matrix(genotypes[, !names(genotypes) %in% c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT")])
  genotypes = as.data.frame(matrix(as.numeric(sub(".*:", "", genotypes)),
                                   nrow = nrow(genotypes), ncol = ncol(genotypes),
                                   dimnames = list(row.names(genotypes), colnames(genotypes))))
  
  return(genotypes)
}

################################################################

args = commandArgs(trailingOnly = TRUE)
RawDir = args[1]
RawFile = args[2]

genotypes = read_genotypes(RawDir, RawFile)
chrs = sub(":.*", "", row.names(genotypes))

for (chr in unique(chrs)){
  genotypes_chr = genotypes[chrs == chr,]
  saveRDS(genotypes_chr, paste(RawDir, paste0(RawFile, "-", chr, ".rds"), sep = "/"))
}
