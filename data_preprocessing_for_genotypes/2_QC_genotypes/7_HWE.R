#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
hwe_threshold = as.numeric(args[2])

hwe = read.table(paste0(prefix, ".hwe"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)
hwe = hwe[hwe[, "TEST"] == "UNAFF",]

#############
# 2) PREPARATION

hwe_threshold_log = -log10(hwe_threshold)
hwe[, "-logP"] = -log10(hwe[, "P"])

observed = -log10(sort(hwe[, "P"]))
expected = -log10(seq_along(hwe[, "P"])/length(hwe[, "P"]+1))

#############
# 3) PLOT

jpeg(paste(dirname(prefix), "plot-HWE-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(hwe[, "-logP"], col = "grey", breaks = 100, xlab = "-log10(HWE exact p-value)", main = "")
abline(v = hwe_threshold_log, col = "red", lty = "dashed")
dev.off()

###

jpeg(paste(dirname(prefix), "plot-HWE-QQ.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(expected, observed, xlab = "-log10(expected p-value)", ylab = "-log10(observed p-value)", main = "")
abline(0, 1, lty = "dashed")
abline(v = hwe_threshold_log, col = "red", lty = "dashed")
dev.off()
