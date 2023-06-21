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
#    path to directory containing fast_bin_perm.r
#    column in phenotype file
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
#trait.forShuf = args[16]
trait.col = args[16:length(args)]
trait.all = traits[trait.col]
if (length(args) == 16) {
        # Convert trait.all into an array
        trait.all = as.matrix(trait.all)
}
#trait.allForShuf = traits[[trait.forShuf]]
valid = c()
for (i in 1:nrow(trait.all)) {
  # Iterate through the rows of the matrix and find those with no NAs
  NAPresent = FALSE
  for (j in 1:ncol(trait.all)) {
    # Iterate through the columns of the matrix and check if each entry is an NA
    if (is.na(trait.all[i,j])) {
      # The entry is an NA
      NAPresent = TRUE
      break
    }
  }
  if (NAPresent == FALSE) {
    # No NAs in current row
    valid = c(valid, i)
  }
}
trait = as.matrix(trait.all[valid,])
if (length(valid) == 0) {
  # No rows with values from all species
  print("Warning: No species with phenotype annotations for all phenotypes.")
}
#traitForShuf = trait.allForShuf[valid]
species.spaces = traits$Species.Name[valid]
trait.species = gsub(" ", "_", species.spaces)
row.names(trait) = trait.species
traitForShuf = trait[,1]

#Read activity predictions
preds = read.csv(file = args[7], header = F, sep = "\t")
names(preds)[1] = "OCR"
te = read.csv(file = args[8], header=F)
pred.species = gsub(" ", "_", te$V1)
if (length(pred.species)+1 != ncol(preds)) {
  print("Warning: Number of species names does not match number of species.")
}
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
  source(paste(args[15], "/fast_bin_perm.r", sep=""))
}

max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * num_shuffles
t = 0.5 * length(common.species)

enh.names = character(n)
p.vals = double(n)
coeffs = matrix(nrow=n, ncol=length(args) - 15)
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
  int.trait = as.matrix(trait[int.species, ])
  int.traitForShuf = traitForShuf[int.species]
  sum.trait = sum(int.traitForShuf)
  if (((l >= t) && (sum.trait > 0)) && (sum.trait < l)) { #50%
    int.preds = good.preds[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      leafMap=makeLeafMap(int.tree.di)
      fg.species = names(int.traitForShuf[which(int.traitFoShuf == 1)])
      bg.species = names(int.traitForShuf[which(int.traitForShuf == 0)])
      fg.leaf.count = length(fg.species)
      fg.internal.count = countInternal(int.tree.di, leafMap, fg.species)
      rate.matrix=ratematrix(int.tree.di, int.traitForShuf)
      int.traitForShuf.real = int.traitForShuf
      names(int.traitForShuf.real) = int.species
    }
    if (length(args) > 16) {
      # Other traits should be used as additional covariates
      for (j in 2:ncol(int.trait)) {
        # Add the other traits as covariates
        int.preds = cbind(as.double(int.preds), as.double(int.trait[,j]))
      }
    } else {
      # Convert int.preds into a double
      int.preds = as.double(int.preds)
    }
    for (f in 1:num_shuffles) {
      if (random) {
        fg.species.shuffled = fastSimBinPhenoVec(int.tree.di, tips=fg.leaf.count, fg.internal.count, rm=rate.matrix, leafBitMaps=leafMap)
        int.trait[,1] = double(l)
        names(int.trait) = int.species
        int.trait[,1][fg.species.shuffled] = 1
      }
      X = int.preds
      Y = int.trait[,1]
      dat <- data.frame(X = X, Y = Y)#, row.names = int.species)
      m <- tryCatch(
        {
        phyloglm(Y ~ X, data = dat, phy=int.tree.di,  method = "logistic_MPLE")
        },
        error=function(e) {
          #message('An Error Occurred')
          #print(e)
          print(name)
          return(NULL)
        })
      if (!is.null(m)) {
        enh.names[index] = name
        m.coeff = summary(m)$coefficients
        p.vals[index] = m.coeff[8 + 3*(length(args) - 16)]
        coeffs[index,] = m.coeff[2:(length(args)-14)]
        index = index + 1
      }
    }
  }
}
proc.time() - ptm
options(warn = 1)

#Output
dat = data.frame(OCR = enh.names[1:index-1], Pvalue = p.vals[1:index-1])
for (i in 1:ncol(int.trait)) {
  # Iterate through the additional coefficients and add them to the data frame
  dat = cbind(dat, Coeff = coeffs[1:index-1,i])
}
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)
