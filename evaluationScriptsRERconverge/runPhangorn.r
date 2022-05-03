args = commandArgs(trailingOnly=TRUE)
inputFileName = args[1]
outputFileName = paste(inputFileName, ".tree.txt", sep="")
fullName = basename(inputFileName)
name = sub("\\.\\w*$", "", fullName)
library("RERconverge")
tree.res = estimatePhangornTree(alnfile = inputFileName, treefile="/ocean/projects/bio200034p/ikaplow/PredictionsToPhenotypes/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree", type="DNA", submodel="GTR", k=4)
fc = file(outputFileName, "wt")
writeLines(paste(name, write.tree(tree.res$tree), sep = "\t"), fc)
