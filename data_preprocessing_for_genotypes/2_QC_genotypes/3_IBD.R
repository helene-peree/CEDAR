#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]

genome_threshold = 0.185
genome = read.table(paste0(prefix, ".genome"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
imiss = read.table(paste0(prefix, ".imiss"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) PREPARATION

replacement = "SYSCID00|SYSCID0|SYSCID|CROHN00|CROHN0|CROHN"

genome_subset = genome[which(genome[, "PI_HAT"] > genome_threshold),]
genome_text = genome[which(genome[, "PI_HAT"] > 0.08),]
genome_text[, "IID1b"] = sub(replacement, "", genome_text[, "IID1"])
genome_text[, "IID2b"] = sub(replacement, "", genome_text[, "IID2"])

#############
# 2) PLOTS

jpeg(paste(dirname(prefix), "plot-IBD.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(genome[, "PI_HAT"], xlab = "Pairs of individuals", ylab = "Identity by descent (IBD)",
     ylim = c(0, 1))
abline(h = genome_threshold, col = "RED", lty = "dashed")
text(row.names(genome_text), genome_text[, "PI_HAT"]-0.025, labels = paste(genome_text[, "IID1b"], genome_text[, "IID2b"], sep = "_"))
dev.off()

###

jpeg(paste(dirname(prefix), "plot-IBD-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(genome[, "PI_HAT"], col = "grey", breaks = 100, xlab = "Identity by descent (IBD)", main = "")
abline(v = genome_threshold, col = "red", lty = "dashed")
dev.off()

#############
# 3) EXPORT EXCLUDED INDIVIDUALS

genome_subset = genome_subset[!(genome_subset[, "IID1"] == "SYSCID036" & genome_subset[, "IID2"] == "CROHN038"),]
genome_subset = genome_subset[!(genome_subset[, "IID1"] == "CROHN038" & genome_subset[, "IID2"] == "SYSCID036"),]

genome_excluded = NULL

for (pair in 1:nrow(genome_subset)){
  F_MISS_ind_1 = imiss[imiss[, "IID"] == genome_subset[pair, "IID1"], "F_MISS"]
  F_MISS_ind_2 = imiss[imiss[, "IID"] == genome_subset[pair, "IID2"], "F_MISS"]
  
  if (F_MISS_ind_1 > F_MISS_ind_2){
    excluded_individual = genome_subset[pair, c("FID1", "IID1")]
    
  } else {
    excluded_individual = genome_subset[pair, c("FID2", "IID2")]
  }
  names(excluded_individual) = c("FID", "IID")
  genome_excluded = rbind(genome_excluded, excluded_individual)
}

genome_excluded = unique(genome_excluded[order(genome_excluded["IID"]),])

write.table(genome_excluded, file = paste(dirname(prefix), "fail-qc-genome.txt", sep = "/"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
