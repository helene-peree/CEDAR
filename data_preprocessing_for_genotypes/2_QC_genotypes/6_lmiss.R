#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
lmiss_threshold = as.numeric(args[2])

lmiss = read.table(paste0(prefix, ".lmiss"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) PLOT

jpeg(paste(dirname(prefix), "plot-lmiss-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(lmiss[, "F_MISS"], col = "grey", breaks = 50, xlab = "Proportion of missing data", main = "")
abline(v = lmiss_threshold, col = "red", lty = "dashed")
dev.off()
