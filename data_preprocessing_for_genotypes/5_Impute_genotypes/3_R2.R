#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]
accuracy_threshold = as.numeric(args[2])

accuracy = read.table(paste0(prefix, ".r2"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) EXCLUSION

accuracy[, "R2"] = as.numeric(sapply(accuracy[, "INFO"], function(x){sub("R2=", "", strsplit(x, ";")[[1]][3])}))
accuracy_excluded = accuracy[accuracy[, "R2"] <= accuracy_threshold, "ID"]

print(length(accuracy_excluded))

#write.table(accuracy_excluded, file = paste(dirname(prefix), "fail-qc-accuracy.txt", sep = "/"),
#            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#############
# 3) PLOT

jpeg(paste(dirname(prefix), "plot-R2-hist.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
hist(accuracy[, "R2"], col = "grey", breaks = 40, xlab = "Imputation accuracy (R2)", main = "")
abline(v = accuracy_threshold, col = "red", lty = "dashed")
dev.off()
