library(ggplot2)

#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
imiss_threshold = as.numeric(args[2])

imiss = read.table(paste0(prefix, ".imiss"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
het = read.table(paste0(prefix, ".het"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) PREPARATION

replacement = "SYSCID00|SYSCID0|SYSCID|CROHN00|CROHN0|CROHN"
replacement2 = "SYSCID00|SYSCID0|SYSCID"

imiss[, "logF_MISS"] = log10(imiss[, "F_MISS"])
imiss[, "Individuals"] = as.integer(rownames(imiss))
imiss_threshold_log = log10(imiss_threshold)

het[, "meanHet"] = (het[, "N(NM)"] - het[, "O(HOM)"])/het[, "N(NM)"]
threshold_het_pos = 0.3165
threshold_het_neg = 0.300
threshold_inbreeding_pos = 0.035
threshold_inbreeding_neg = -0.017

imiss_subset = imiss[imiss[, "F_MISS"] > imiss_threshold,]
het_subset = het[het[, "meanHet"] < threshold_het_neg | het[, "meanHet"] > threshold_het_pos,]
imiss_subset_het = imiss[imiss[, "IID"] %in% het_subset[, "IID"] & imiss[, "F_MISS"] < imiss_threshold,]
het_subset_imiss = het_subset[het_subset[, "IID"] %in% imiss_subset_het[, "IID"],]
het_subset_inbreeding = het[het[, "F"] < threshold_inbreeding_neg | het[, "F"] > threshold_inbreeding_pos,]

colors = densCols(imiss[, "logF_MISS"], het[, "meanHet"])

#############
# 3) PLOTS

jpeg(paste(dirname(prefix), "plot-imiss.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(imiss[, "F_MISS"], xlab = "Individuals", ylab = "Missing call rate", main = "",
     col = ifelse(!grepl("old", imiss[, "IID"]), "black", "blue"))
abline(h = imiss_threshold, col = "red", lty = "dashed")
abline(h = 0.02, col = "red", lty = "dashed")
text(row.names(imiss_subset), imiss_subset[, "F_MISS"]-0.005, labels = sub(replacement, "", imiss_subset[, "IID"]),
     col = ifelse(!grepl("old", imiss_subset[, "IID"]), "black", "blue"))
dev.off()

###

jpeg(paste(dirname(prefix), "plot-imiss-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(imiss[, "F_MISS"], col = "grey", breaks = 20, xlab = "Missing call rate", main = "")
abline(v = imiss_threshold, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-imiss-het.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(imiss[, "F_MISS"], het[, "meanHet"],
     xlab = "Missing call rate", ylab = "Heterozygosity rate",
     col = colors, axes = F, ylim = c(0.2, 0.5), pch = 20)
axis(1, at = seq(0, 0.36, by = 0.02))
axis(2, at = seq(0, 0.5, by = 0.05))
abline(v = imiss_threshold, col = "red", lty = "dashed")
abline(h = threshold_het_pos, col = "red", lty = "dashed")
abline(h = threshold_het_neg, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-imiss-het-log.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(imiss[, "logF_MISS"], het[, "meanHet"],
     xlab = "Proportion of missing data (log10 scale)", ylab = "Heterozygosity rate",
     col = colors, axes = F, ylim = c(0.2, 0.5), pch = 20) #xlim = c(-3, 0)
axis(1, at = c(-3, -2, imiss_threshold_log, -1, 0), labels = c(0.001, 0.01, imiss_threshold, 0.1, 1))
axis(2, at = seq(0, 0.5, by = 0.05))
abline(v = imiss_threshold_log, col = "red", lty = "dashed")
abline(h = threshold_het_pos, col = "red", lty = "dashed")
abline(h = threshold_het_neg, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-imiss-het-log-zoom.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(imiss[, "logF_MISS"], het[, "meanHet"],
     xlab = "Proportion of missing data (log10 scale)", ylab = "Heterozygosity rate",
     xaxt = "n", xlim = c(-3, imiss_threshold_log),
     ylim = c(min(het[, "meanHet"]-0.0035), max(het[, "meanHet"])))
axis(1, at = c(-3, -2, imiss_threshold_log), labels = c(0.001, 0.01, imiss_threshold))
text(imiss_subset_het[, "logF_MISS"], het_subset_imiss[, "meanHet"]-0.003, sub(replacement2, "", het_subset_imiss[, "IID"]))
abline(h = threshold_het_pos, col = "red", lty = "dashed")
abline(h = threshold_het_neg, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-het.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(het[, "meanHet"], xlab = "Individuals", ylab = "Heterozygosity rate",
     xlim = c(min(as.numeric(row.names(het))), max(as.numeric(row.names(het)))+20),
     ylim = c(min(het[, "meanHet"]-0.0035), max(het[, "meanHet"])))
text(row.names(het_subset), het_subset[, "meanHet"]-0.003, sub(replacement2, "", het_subset[, "IID"]))
abline(h = threshold_het_pos, col = "red", lty = "dashed")
abline(h = threshold_het_neg, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-inbreeding.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(het[, "F"], xlab = "Individuals", ylab = "F coefficient", main = "",
     xlim = c(min(as.numeric(row.names(het))), max(as.numeric(row.names(het)))+20),
     ylim = c(min(het[, "F"]-0.015), max(het[, "F"])))
text(row.names(het_subset_inbreeding), het_subset_inbreeding[, "F"]-0.01, sub(replacement2, "", het_subset_inbreeding[, "IID"]))
abline(h = threshold_inbreeding_pos, col = "red", lty = "dashed")
abline(h = threshold_inbreeding_neg, col = "red", lty = "dashed")
dev.off()

#############
# 4) EXPORT EXCLUDED INDIVIDUALS

imiss_excluded = imiss[imiss[, "F_MISS"] >= imiss_threshold, c("FID", "IID")] 

write.table(imiss_excluded, file = paste(dirname(prefix), "fail-qc-imiss.txt", sep = "/"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
