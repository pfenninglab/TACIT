#Script to preform phyloglm analysis on OCR activity predictions. This script is for command-line use.
#Daniel Schaffer, Pfenning Lab
#Takes several required positional arguments, in order: 
#    tree file (Newick), 
#    OCR activity predictions matrix file (no header, first column is OCR names, tab-seperated), 
#    file listing species corresponding to predictions matrix columns (one per line, in the format Genus species), 
#    phenotype file (CSV, species names in a column called "Species Name" in the format Genus species)
#    output file name "template" (see below), 
#    I, output file number and inital line in predictions matrix (see below), 
#    J, step size in predictions matrix (see below), 
#    K, number of permulations per OCR (0 for true data / no permulations), 
#    S, random seed to use, 
#    column in phenotype file
#    path to directory containing fast_bin_perm.r (not required if K=0)
#E.g. if putput file "template" given is /path/to/foo.csv, output will be in (I,S as above) /path/to/foo_rI_S.csv
#Template must end in .csv
#Will apply phylolm to K permulations for OCRs on lines I, I+J, I+2J, ... until end of matrix is reached

library(ape)     #Phylogenetic tree processing
library(phylolm) #Phylogeny-corrected correlation


args <- commandArgs()
seed = as.integer(args[14])
set.seed(seed)

tree <- read.tree(file = args[6]) #Change to read.nexus for a nexus-format tree

#Read phenotype data
traits = read.csv(file= args[9])
trait.col = args[15]
trait.all = traits[[trait.col]]
valid = which(trait.all %in% c(0,1))
trait = as.numeric(as.character(trait.all[valid]))
species.spaces = traits$Species.Name[valid]
trait.species = gsub(" ", "_", species.spaces)
names(trait) = trait.species

#Read activity predictions
preds = read.csv(file = args[7], header = F, sep = "\t")
names(preds)[1] = "OCR"
te = read.csv(file = args[8], header=F)
pred.species = gsub(" ", "_", te$V1)
names(preds)[2:(length(pred.species)+1)] = pred.species

common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
te = which(trait.species %in% common.species)
tree.common = keep.tip(tree, common.species)


########

#Setup
row_init = as.integer(args[11])
row_step = as.integer(args[12])
num_shuffles = as.integer(args[13])

random = T
if (num_shuffles == 0) {
  random=F
  num_shuffles = 1
} else {
  source(paste(args[16], "/fast_bin_perm.r", sep=""))
}

max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * num_shuffles
t = 0.5 * length(common.species)

enh.names = character(n)
p.vals = double(n)
coeffs = double(n)
index = 1

#Iterate & run phyloglm
options(warn = -1) #suppress warning from phyloglm that boundaries of parameters are reached
ptm <- proc.time()
for (i in 0:max_iter) {
  e = row_init + i*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  l = length(int.species)
  int.trait = trait[int.species]
  sum.trait = sum(int.trait)
  if (l >= t && sum.trait > 0 && sum.trait < l) { #50%
    int.preds = good.preds[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      leafMap=makeLeafMap(int.tree.di)
      fg.species = names(int.trait[which(int.trait == 1)])
      bg.species = names(int.trait[which(int.trait == 0)])
      fg.leaf.count = length(fg.species)
      fg.internal.count = countInternal(int.tree.di, leafMap, fg.species)
      rate.matrix=ratematrix(int.tree.di, int.trait)
    }
    for (f in 1:num_shuffles) {
      if (random) {
        fg.species.shuffled = fastSimBinPhenoVec(int.tree.di, tips=fg.leaf.count, fg.internal.count, rm=rate.matrix, leafBitMaps=leafMap)
        int.trait = double(l)
        names(int.trait) = int.species
        int.trait[fg.species.shuffled] = 1
      }
      dat <- data.frame(Y=int.trait, X=as.double(int.preds), row.names = int.species)
      m <- phyloglm(Y ~ X, data = dat, phy=int.tree.di,  method = "logistic_MPLE")
      enh.names[index] = name
      m.coeff = summary(m)$coefficients
      p.vals[index] = m.coeff[8]
      coeffs[index] = m.coeff[2]
      index = index + 1
    }
  }
}
proc.time() - ptm
options(warn = 1)

#Output
dat = data.frame(OCR = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)


