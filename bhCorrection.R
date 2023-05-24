#Used to add a column 'bh' to adjust p value with bh
#Usage1:
# Rscript bhCorrection.R input.csv output.csv
#Usage2:
# Rscript bhCorrection.R input.csv output.csv Pvalue

args <- commandArgs(trailingOnly = TRUE)
input <- args[1]
output <- args[2]
# optional adjust column name
if (length(args) >= 3) {
  pName <- args[3]
} else {
  pName <- "Exp_Pvalue"
}
data <- read.csv(input)
pval <- data[[pName]]
bh = p.adjust(pval, method = "BH")
data$bh=bh
# sort table by adjusted P
sorted <- data[order(data$bh),]
write.table(sorted,output,sep=",",row.names = FALSE,quote = FALSE)