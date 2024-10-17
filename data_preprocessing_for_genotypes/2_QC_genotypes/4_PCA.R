#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]

hapmap3 = read.table(sub("Output.*", "HapMap_article/relationships_w_pops_121708.txt", prefix), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
eigenvec = read.table(paste0(prefix, ".eigenvec"), stringsAsFactors = FALSE, check.names = FALSE, header = FALSE,
                      col.names = c("FID", "IID", "PC1", "PC2", "PC3"))

#############
# 2) PREPARATION

hapmap3 = hapmap3[grep("CEU|CHB|JPT|YRI", hapmap3[, "population"]),]
eigenvec[, "population"] = hapmap3[match(eigenvec[, "IID"], hapmap3[, "IID"]), "population"]
eigenvec[is.na(eigenvec[, "population"]), "population"] = "CEDAR_II"

threshold_PC1 = 13
threshold_PC1_pos = round(mean(eigenvec[grep("CEU", eigenvec[, "population"]), "PC1"]) + threshold_PC1*(sd(eigenvec[grep("CEU", eigenvec[, "population"]), "PC1"])), 3)
threshold_PC1_neg = round(mean(eigenvec[grep("CEU", eigenvec[, "population"]), "PC1"]) - threshold_PC1*(sd(eigenvec[grep("CEU", eigenvec[, "population"]), "PC1"])), 3)

eigenvec_subset2 = eigenvec[(eigenvec[, "PC1"] > threshold_PC1_pos | eigenvec[, "PC1"] < threshold_PC1_neg) & eigenvec[, "population"] == "CEDAR_II",]

#############
# 3) PLOT

jpeg(paste(dirname(prefix), "plot-PCA.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(eigenvec[, "PC1"], eigenvec[, "PC2"], col = as.factor(eigenvec[, "population"]),
     pch = 20, xlab = "PC1", ylab = "PC2")
abline(v = threshold_PC1_pos, lty = "dashed")
abline(v = threshold_PC1_neg, lty = "dashed")
text(eigenvec_subset2[, "PC1"], eigenvec_subset2[, "PC2"], eigenvec_subset2[, "IID"])
legend("topleft",legend = unique(as.factor(eigenvec[, "population"])), pch = 20,
       col = unique(as.factor(eigenvec[, "population"])))
dev.off()

#############
# 4) EXPORT EXCLUDED INDIVIDUALS

eigenvec_excluded = eigenvec_subset2[, c("FID", "IID")]

write.table(eigenvec_excluded, file = paste(dirname(prefix), "fail-qc-pca.txt", sep = "/"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
