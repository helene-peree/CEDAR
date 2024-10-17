#############
# 1) IMPORT

args = commandArgs(trailingOnly=TRUE)
prefix = args[1]

sex = read.table(paste0(prefix, ".sexcheck"), stringsAsFactors = FALSE, check.names = FALSE, header = TRUE)

#############
# 2) PREPARATION

replacement = "SYSCID00|SYSCID0|SYSCID|CROHN00|CROHN0|CROHN"
sex_subset = sex[(sex[, "F"] > 0.2 & sex[, "F"] < 0.8),]

#############
# 3) PLOT

jpeg(paste(dirname(prefix), "plot-sex.jpg", sep = "/"), type = "cairo", res = 300, units = "in", width = 7, height = 7)
plot(sex[, "F"], xlab = "Individuals", ylab = "F coefficient", main = "",
     pch = ifelse(sex[, "STATUS"] == "OK", 1, 19),
     col = ifelse(sex[, "STATUS"] == "OK", "black", ifelse(sex[, "PEDSEX"] == 0, "black", "red")))
abline(h = 0.2, col = "red", lty = "dashed")
abline(h = 0.8, col = "red", lty = "dashed")
text(10, 0.15, "females")
text(5, 0.85, "males")
if (nrow(sex_subset) != 0){
  text(row.names(sex_subset), sex_subset[, "F"]-0.025, labels = sub(replacement, "", sex_subset[, "IID"]))
}
dev.off()

#############
# 3) EXPORT EXCLUDED INDIVIDUALS

sex_excluded = sex[sex[, "STATUS"] == "PROBLEM", c("FID", "IID")] 

write.table(sex_excluded, file = paste(dirname(prefix), "fail-qc-sexcheck.txt", sep = "/"),
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
