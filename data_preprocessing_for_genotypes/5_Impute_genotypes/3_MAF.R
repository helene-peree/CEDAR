#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
maf_threshold = as.numeric(args[2])

maf = read.table(paste0(prefix, ".frq"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) EXCLUSION

maf_excluded = maf[maf[, "MAF"] < maf_threshold | maf[, "MAF"] > (1-maf_threshold), "SNP"]

print(length(maf_excluded))

write.table(maf_excluded, file = paste(dirname(prefix), "fail-qc-maf.txt", sep = "/"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#############
# 3) PLOT

jpeg(paste(dirname(prefix), "plot-MAF-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(maf[, "MAF"], col = "grey", breaks = 40, xlab = "Minor allele frequency (MAF)", main = "")
abline(v = maf_threshold, col = "red", lty = "dashed")
abline(v = 1-maf_threshold, col = "red", lty = "dashed")
dev.off()
