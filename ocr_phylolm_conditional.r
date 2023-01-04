#Script to preform phylolm analysis on OCR activity predictions. This script is for command-line use.
#Daniel Schaffer, Pfenning Lab
#Takes several required positional arguments, in order: 
#    tree file (Newick), 
#    OCR activity predictions matrix file (no header, first column is OCR names, tab-seperated), 
#    file listing species corresponding to predictions matrix columns (one per line, in the format Genus species), 
#    phenotype file (CSV, species names in a column called "species binomial" in the format genus_species)
#    output file name "template" (see below), 
#    I, output file number and inital line in predictions matrix (see below), 
#    J, step size in predictions matrix (see below), 
#    (changed) The name of a CSV with columns Enhancer and Missing_Trials
#         that specifies how many permulations (K) to do for each OCR
#    S, random seed to use, 
#    column in phenotype file
#E.g. if putput file "template" given is /path/to/foo.csv, output will be in (I,S as above) /path/to/foo_rI_S.csv
#Template must end in .csv
#Will apply phylolm to K permulations for OCRs on lines I, I+J, I+2J, ... until end of matrix is reached

library(ape)     #Phylogenetic tree processing
library(geiger)  #Brownian motion simulation
library(stringr) #Species name manipulation
library(phylolm) #Phylogeny-corrected correlation


#Slightly modified from RERconverge package
#'Generates a permulated continuous phenotype given an observed continuous phenotype and a phylogeny
#' @param namedvec A named numeric vector with phenotype values for each speices
#' @param treewithbranchlengths A rooted phylogenetic tree with the same species as namedvec and branch lengths representing average evolutionary rate.
#' @param rm A precomputed rate matrix
#' @return A vector with permulated phenotype values
#' @export
simpermvec=function (namedvec, treewithbranchlengths, rm=NULL) 
{
  if(is.null(rm)) {
    rm = ratematrix(treewithbranchlengths, namedvec)
  }
  sims = sim.char(treewithbranchlengths, rm, nsim = 1)
  nam = rownames(sims)
  s = as.data.frame(sims)
  simulatedvec = s[, 1]
  names(simulatedvec) = nam
  simsorted = sort(simulatedvec)
  realsorted = sort(namedvec)
  l = length(simsorted)
  c = 1
  while (c <= l) {
    simsorted[c] = realsorted[c]
    c = c + 1
  }
  simsorted
}

args = commandArgs()
seed = as.integer(args[14])
set.seed(seed)

tree = read.tree(file = args[6]) #Change to read.nexus for a nexus-format tree

#Read phenotype data
traits = read.csv(file= args[9])
trait.col = args[15]
trait.all = traits[[trait.col]]
valid = which(trait.all != "NA") 
trait = trait.all[valid]
species.lower = traits$species.binomial[valid]
trait.species = str_to_title(species.lower)
names(trait) = trait.species

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
enh_details = read.csv(file = args[13], header = T)
enh_shuffles = enh_details$Missing_Trials
names(enh_shuffles) = enh_details$Enhancer
enh_coeffs = enh_details$Coeff
names(enh_coeffs) = enh_details$Enhancer

max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * max(enh_shuffles)

enh.names = character(n)
p.vals = double(n)
coeffs = double(n)
index = 1

#Iterate & run phylolm
options(warn = -1) #suppress warning from phylolm
ptm <- proc.time()
for (i in 0:max_iter) {
  e = row_init + i*row_step
  name = as.character(preds[e, 1])
  if (name %in% names(enh_shuffles)) {
    num_shuffles = enh_shuffles[name]
    orig_coeff_sign = sign(enh_coeffs[name])  

    cur.preds = preds[e, 2:(length(pred.species)+1)]
    good.preds = cur.preds[which(cur.preds != -1)]
    int.species = intersect(names(good.preds), common.species)
    int.preds = good.preds[int.species]
    int.trait = trait[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    rate.matrix=ratematrix(int.tree.di, int.trait)
    int.trait.real = int.trait
    names(int.trait.real) = int.species
  
    for (f in 1:num_shuffles){ 
	  repeat {
        int.trait = simpermvec(int.trait.real, int.tree.di, rm=rate.matrix)
        dat <- data.frame(X = as.double(int.preds), Y = int.trait)
        m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")
        m.coeff = summary(m)$coefficients
	    if (sign(m.coeff[2]) == orig_coeff_sign) {
	      enh.names[index] = name
          p.vals[index] = m.coeff[8]
          coeffs[index] = m.coeff[2]
          index = index + 1
		  break
	    }
      }
    }
  }	
}
proc.time() - ptm
options(warn = 1)

#Output
dat = data.frame(OCR = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)
