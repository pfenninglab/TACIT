#Script to preform PGLMM analysis on enhancer activity predictions. This script is for Windows/RStudio or command-line use, depedning on which lines are (un)commented.

library(ape)
library(geiger)
library(stringr)
library(phylolm)

#Copied from RERconverge
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


#Positional Arguments: 
#    tree, 
#    predictions, 
#    species of prediction columns, 
#    traits, 
#    output file template, 
#    output file number and inital line in preds file, 
#    step size in preds file, 
#    number of trials, 
#    seed to use, 
#    trait column
#E.g. if putput file given is /path/to/foo.csv and number is i, output is in /path/to/foo_ri.csv
#Output file # must be the same as the line to start in the input file (laziness)
#Template must end in .csv
#What trait to use is currently still hardcoded
# given last three arguments i,j,k,
#does k shiffles for enhancer i, i+j, i+2j, ... until end of file is reached


args <- commandArgs()
seed = as.integer(args[14])
set.seed(seed)
#print(args)
#print(seed)

#ATTN: Change method to read.tree for newick tree and read.nexus for nexus tree
tree <- read.tree(file = args[6])
#tree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")

##For reading panthiera
#traits = read.csv(file= args[9])
#traits = read.csv("panthiera.csv")
#trait.col = args[15]
#trait.col = "X9.1.GestationLen.d"
#trait.all = traits[[trait.col]]
#iv = which(trait.all > 0)
#trait = trait.all[iv]
#species.spaces = traits$Species[iv]
#trait.species = sub(" ", "_", species.spaces) #Alert! duplicates. Not sure if matters.

##For reading zoonomia sheet
traits = read.csv(file= args[9])
#traits = read.csv("teeling_longevity_2-12-21.csv")
#traits = read.csv("zoonomia_phenotypes_4-16-21.csv")
trait.col = args[15]
#trait.col = "LQ"
#trait.col= "Brain.resid"
#trait.col = "Resting.body.temperature"
trait.all = traits[[trait.col]]
#iv = which(trait.all > 0) 
iv = which(trait.all != "NA") 
trait = trait.all[iv]
species.lower = traits$species.binomial[iv]
trait.species = str_to_title(species.lower)
names(trait) = trait.species

clades.all = traits$Taxonomic.lineage
names(clades.all) = str_to_title(traits$species.binomial)


#Possible TODO - replace reading with reading one line at a time. Might reduce memory usage by each run
preds = read.csv(file = args[7], header = F, sep = "\t")
#preds = read.csv(file = "liver\\MoRaMaCoPiEnhancersNRSummit.txt", header = F, sep = "\t")
#preds = read.csv(file = "cortex\\MoMaRaBaEnhancersNRSummit.txt", header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = args[8], header=F)
#te = read.csv(file = "BoreoeutheriaSpeciesNamesNew.txt", header=F)
pred.species = gsub(" ", "_", te$V1)
names(preds)[2:(length(pred.species)+1)] = pred.species

common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
te = which(trait.species %in% common.species)
tree.common = keep.tip(tree, common.species)


########

row_init = as.integer(args[11])
row_step = as.integer(args[12])
num_shuffles = as.integer(args[13])
#row_init = 1
#row_step = 1
#num_shuffles = 0 #Should be 0 for non-random

random = T
if (num_shuffles == 0) {
  random=F
  num_shuffles = 1
}


max_iter = (nrow(preds)-row_init) %/% row_step
n = (max_iter + 1) * num_shuffles
t = 0.5 * length(common.species)

enh.names = character(n)
p.vals = double(n)
coeffs = double(n)
index = 1

options(warn = -1) #suppress warning from phylolm
ptm <- proc.time()
for (m in 0:max_iter) {
  e = row_init + m*row_step
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  l = length(int.species)
  if (l >= t) { #50%
    int.preds = good.preds[int.species]
    int.trait = trait[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
    if (random) {
      rate.matrix=ratematrix(int.tree.di, int.trait)
      int.trait.real = int.trait
      names(int.trait.real) = int.species
    }
  
    for (f in 1:num_shuffles) {
      print("hi")
      if (random) {
        int.trait = simpermvec(int.trait.real, int.tree.di, rm=rate.matrix)
      }
      dat <- data.frame(X = as.double(int.preds), Y = int.trait)
      m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")
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

int.clades = clades.all[rownames(dat)]
color = rep("gray", nrow(dat))
color[which(grepl("Laurasiatheria", int.clades))] = "red"
color[which(grepl("Glires", int.clades))] = "blue"
color[which(grepl("Primates", int.clades))] = "darkgreen"

svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,".svg",sep=""), width=8, height=6)

plot(dat$X, dat$Y, col=color, pch=19,  xlab = "Enhancer activity prediction", ylab = "Brain Size", main=name, cex.lab=1.5, cex.main = 2, cex.axis=1.5)


abline(a=m$coefficients[1], b=m$coefficients[2], lwd=2)
dev.off()

color.scale = c("#164d98", "#3f4794", "#58408e", "#6b3986", "#7b307b", "#87276e", "#901f60", "#961951", "#991741", "#991b32")
color = rep("gray", nrow(dat))
for (i in 1:nrow(dat)) {
  sp = int.tree$tip.label[i]
  pred = dat[sp,1]
  color[i] = color.scale[as.integer(pred*10) + 1]
}

plot(int.tree, tip.color=color)



X_range <- seq(from=min(dat$X), to=max(dat$X), by=.01)
y_logits <- m$coefficients[1] + m$coefficients[2]*X_range
y_probs <- exp(y_logits)/(1 + exp(y_logits))
points(X_range, y_probs, ylim=c(0,1), type="l", lwd=2, lty=1)
dev.off()



for (i in 1:35) {
  print(i)
  e = sigrows[i]
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
    int.preds = good.preds[int.species]
    int.trait = trait[int.species]
    int.tree = keep.tip(tree.common, int.species)
    int.tree.di = multi2di(int.tree)
      dat <- data.frame(X = as.double(int.preds), Y = int.trait)

      m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")
      
      int.clades = clades.all[rownames(dat)]
      color = rep("gray", nrow(dat))
      laur = which(grepl("Laurasiatheria", int.clades))
      glir = which(grepl("Glires", int.clades))
      euar = which(grepl("Primates", int.clades))
      color[laur] = "#ff8000"
      color[glir] = "#0000ff"
      color[euar] = "#008000"
      
      shape=rep(19, nrow(dat))
      ape = which(grepl("Hominoidea", int.clades))
      #cet = which(grepl("Cetacea", int.clades))
      cet = which(grepl("Chiroptera", int.clades))
      shape[ape] = 15
      shape[cet] = 17
      
      laur.line = phylolm(Y ~ X, data = dat[laur,], phy=int.tree.di, model = "BM")
      if (length(glir) > 0) {
      glir.line = phylolm(Y ~ X, data = dat[glir,], phy=int.tree.di, model = "BM")
      }
      euar.line = phylolm(Y ~ X, data = dat[euar,], phy=int.tree.di, model = "BM")
      
      
      #svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_v3.svg",sep=""), width=8, height=8)
      svg(file = paste("liver\\lq\\liver_lq_",name,"_v3.svg",sep=""), width=8, height=8)
      #change ylim depending on trait - for brain size, -0.75 - 0.75
      #plot(dat$X, dat$Y, col=color, pch=shape, xlim=c(0, 1), ylim=c(-0.75,0.75), xlab="OCR activity prediction", cex=2, ylab = "Brain size residual", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
      plot(dat$X, dat$Y, col=color, pch=shape, xlim=c(0, 1), ylim=c(0,10), xlab="OCR activity prediction", cex=2, ylab = "Longevity Quotient", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
      
      abline(a=m$coefficients[1], b=m$coefficients[2], lwd=4)
      
      abline(a=laur.line$coefficients[1], b=laur.line$coefficients[2], lty=2, col="#ff8000", lwd=2)
      if (length(glir) > 0) {
      abline(a=glir.line$coefficients[1], b=glir.line$coefficients[2], lty=2, col="#0000ff", lwd=2)
      }
      abline(a=euar.line$coefficients[1], b=euar.line$coefficients[2], lty=2, col="#008000", lwd=2)
      
      
      #legend("bottom", x.intersp=1, xjust=0, yjust=0, bty="n", text.width=c(0.14,0.05,0.12), legend=c("Laurasiatheria", "Glires", "Euarchonta"), col=c("red", "blue", "darkgreen"), pch=c(19,19,19), text.col="black", horiz=T)
      #rect(0.2,-0.9,0.8,-0.65)
      
      dev.off()
      
      
}


pv.table = read.csv("PVenhPreds_Euarchontoglires_brainResidSig.csv")
#pv.enh = "mm10_chr13.114793237.114793737"
#pv.enh = "mm10_chr13.114757413.114757913"
#pv.enh = "mm10_chr2.94470700.94471200"
#pv.enh = "mm10_chr1.95762160.95762660"
name = pv.enh
pv.table.trim = pv.table[which(pv.table[[pv.enh]] > -2),]

int.species = pv.table.trim$species.binomial
int.preds = pv.table.trim[[pv.enh]]
names(int.preds) = int.species
int.trait = pv.table.trim$Brain.resid
names(int.trait) = int.species
int.tree = keep.tip(tree.common, int.species)
int.tree.di = multi2di(int.tree)
dat <- data.frame(X = as.double(int.preds), Y = int.trait)
m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")


int.clades = pv.table.trim$Taxonomic.lineage
color = rep("gray", nrow(dat))
laur = which(grepl("Laurasiatheria", int.clades))
glir = which(grepl("Glires", int.clades))
euar = which(grepl("Primates", int.clades))
color[laur] = "#ff8000"
color[glir] = "#0000ff"
color[euar] = "#008000"

shape=rep(19, nrow(dat))
ape = which(grepl("Hominoidea", int.clades))
cet = which(grepl("Cetacea", int.clades))
shape[ape] = 15
shape[cet] = 17

#laur.line = phylolm(Y ~ X, data = dat[laur,], phy=int.tree.di, model = "BM")
if (length(glir) > 0) {
  glir.line = phylolm(Y ~ X, data = dat[glir,], phy=int.tree.di, model = "BM")
}
euar.line = phylolm(Y ~ X, data = dat[euar,], phy=int.tree.di, model = "BM")


svg(file = paste("cortex\\brain_size\\pv_brain-size_",name,"_v3.svg",sep=""), width=8, height=8)

#change ylim depending on trait - for brain size, -0.75 - 0.75
plot(dat$X, dat$Y, col=color, pch=shape, xlim=c(0, 1), ylim=c(-0.75,0.75), xlab="OCR activity prediction", cex=2, ylab = "Brain size residual", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)

abline(a=m$coefficients[1], b=m$coefficients[2], lwd=4)

#abline(a=laur.line$coefficients[1], b=laur.line$coefficients[2], lty=2, col="#ff8000", lwd=2)
if (length(glir) > 0) {
  abline(a=glir.line$coefficients[1], b=glir.line$coefficients[2], lty=2, col="#0000ff", lwd=2)
}
abline(a=euar.line$coefficients[1], b=euar.line$coefficients[2], lty=2, col="#008000", lwd=2)


dev.off()




library("viridis")




#NEW-SYTLE FIGURE
  e = 75902 #2434, 45190, 75902, 89187
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  int.preds = good.preds[int.species]
  int.trait = trait[int.species]
  int.tree = keep.tip(tree.common, int.species)
  int.tree.di = multi2di(int.tree)
  dat <- data.frame(X = as.double(int.preds), Y = int.trait)
  
  m <- phylolm(Y ~ X, data = dat, phy=int.tree.di, model = "BM")
  
  #             5,              2            28,           6                2           1        10              11          20           3             5               11              3                14           10               11             15                  5
  cladeList = c("Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Tylopoda", "Suina", "Ruminantia", "Cetacea", "Chiroptera", "Lagomorpha", "Sciuromorpha", "Hystricomorpha", "Castorimorpha", "Myomorpha", "Strepsirrhini", "Platyrrhini", "Cercopithecoidea", "Hominoidea")
  cladeNums = c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 12, 13, 14, 15)
  #Problems: Tylopoda (camels) - 2, Suina (pig) - 1 - for now, seperate
  
  int.clades = clades.all[rownames(dat)]
  cladeX = rep(0, nrow(dat))
  for (i in 1:nrow(dat)) {
    cladeNum = cladeNums[str_detect(int.clades[i], cladeList)] #uses stringr
    if (length(cladeNum) != 1) {
      print(i)
    }
    cladeX[i] = cladeNum
  }

  color.scale = plasma(49)[1:39]
  
  ii <- cut(dat$Y, breaks = seq(-0.75, 0.75, len = 40), include.lowest = TRUE)
  colors = color.scale[ii]
  

  svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_newv2.svg",sep=""), width=16, height=8)
  set.seed(0)
  plot(jitter(cladeX), dat$X, col=colors, pch=19, xlim=c(1, 15), ylim=c(0,1), xlab="Clade", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)


  
  #legend("bottom", x.intersp=1, xjust=0, yjust=0, bty="n", text.width=c(0.14,0.05,0.12), legend=c("Laurasiatheria", "Glires", "Euarchonta"), col=c("red", "blue", "darkgreen"), pch=c(19,19,19), text.col="black", horiz=T)
  #rect(0.2,-0.9,0.8,-0.65)
  
  dev.off()

  
  
  library(s2dverification)
  
  svg(file = paste("cortex\\brain_size\\plasma_color_bar.svg",sep=""), width=8, height=1)
  
  ColorBar(brks=40, cols=color.scale, vertical=F, bar_limits=c(-0.75,0.75), subsampleg=6)

  dev.off()

  
  whales = grep("Cetacea", int.clades)
  
  svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_whales.svg",sep=""), width=6, height=6)
  plot(dat$Y[whales], dat$X[whales], col=colors[whales], pch=19, ylim=c(0.25, 0.8), xlim=c(-0.75,0.75),xlab="Bran size residual", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
  line = lm(X ~ Y, dat[whales,])
  abline(a=line$coefficients[1], b=line$coefficients[2], lwd=4)
  dev.off()

  whales = grep("Cetacea", int.clades)
  
  svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_whales.svg",sep=""), width=6, height=6)
  plot(dat$Y[whales], dat$X[whales], col=colors[whales], pch=19, ylim=c(0.25, 0.8), xlim=c(-0.75,0.75),xlab="Bran size residual", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
  line = lm(X ~ Y, dat[whales,])
  abline(a=line$coefficients[1], b=line$coefficients[2], lwd=4)
  dev.off()
  


#For solitary

  soltraits = read.csv("cortex\\solitary\\SocialPhenotypes.csv")
  iv = which(soltraits$solitary_living_euarchontoglires != "NA")
  soltrait = soltraits$solitary_living_euarchontoglires[iv]
  names(soltrait) = soltraits$Species.Name[iv]
  common.species = intersect(intersect(pred.species, tree$tip.label), names(soltrait))
  tree.common = keep.tip(tree, common.species)
  #NEW-SYTLE FIGURE
  e = 68574 
  name = as.character(preds[e, 1])
  cur.preds = preds[e, 2:(length(pred.species)+1)]
  
  good.preds = cur.preds[which(cur.preds != -1)]
  int.species = intersect(names(good.preds), common.species)
  int.preds = good.preds[int.species]
  int.trait = soltrait[int.species]
  int.tree = keep.tip(tree.common, int.species)
  int.tree.di = multi2di(int.tree)
  dat <- data.frame(X = as.double(int.preds), Y = int.trait)
  

  #             5,              2            28,           6                2           1        10              11          20           3             5               11              3                14           10               11             15                  5
  cladeList = c("Placeholder", "Placeholder", "Placeholder", "Placeholder", "Placeholder", "Placeholder", "Placeholder", "Lagomorpha", "Sciuromorpha", "Hystricomorpha", "Castorimorpha", "Myomorpha", "Strepsirrhini", "Platyrrhini", "Cercopithecoidea", "Hominoidea", "Scandentia")
  cladeNums = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 10, 11, 12, 13, 14, 15)
  #Problems: Tylopoda (camels) - 2, Suina (pig) - 1 - for now, seperate
  
  int.clades = clades.all[rownames(dat)]
  
  names(int.clades)[9] = "Capromys_pilorides"
  int.clades[9] = "Hystricomorpha"
  int.clades[24] = "Hystricomorpha"
  names(int.clades)[24] = "Cuniculus_paca"
  names(int.clades)[77] = "Tupaia_tana"
  int.clades[77] = "Scandentia"
  
  #pv
  names(int.clades)[44] = "Tupaia_tana"
  int.clades[44] = "Scandentia"
  
  
  cladeX = rep(0, nrow(dat))
  for (i in 1:nrow(dat)) {
    cladeNum = cladeNums[str_detect(int.clades[i], cladeList)] #uses stringr
    if (length(cladeNum) != 1) {
      print(i)
    }
    cladeX[i] = cladeNum
  }
  
  color.scale = plasma(49)[1:39]
  
  col0 = color.scale[7]
  col1 = color.scale[39-6]
  colors = rep(col0, nrow(dat))
  colors[which(dat$Y == 1)] = col1

  
  svg(file = paste("cortex\\solitary\\cortex_solitary_",name,"_new.svg",sep=""), width=16, height=8)
  set.seed(0)
  plot(jitter(cladeX), dat$X, col=colors, pch=19, xlim=c(1, 15), ylim=c(0,1), xlab="Clade", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
  
  
  
  #legend("bottom", x.intersp=1, xjust=0, yjust=0, bty="n", text.width=c(0.14,0.05,0.12), legend=c("Laurasiatheria", "Glires", "Euarchonta"), col=c("red", "blue", "darkgreen"), pch=c(19,19,19), text.col="black", horiz=T)
  #rect(0.2,-0.9,0.8,-0.65)
  
  dev.off()
  
  
#PV stuff for above
  pv.preds = read.csv(file = "cortex\\solitary\\PV_MoHu_AvgPredictionsMultispeciesPVmodel_euarchontogliresOnly_OrthIn50percentOfSpecies_noHeader.txt", header = F, sep = "\t")
  names(pv.preds)[1] = "Enhancer"
  te = read.csv(file = "cortex\\solitary\\speciesNames_PV_euarchontoglires.txt", header=F)
  pv.pred.species = te$V1
  names(pv.preds)[2:(length(pv.pred.species)+1)] = pv.pred.species
  
  common.species = intersect(intersect(pv.pred.species, tree$tip.label), names(soltrait))
  tree.common = keep.tip(tree, common.species)
  name = "mm10_chr5.134485808.134486308"
  e = 25206
  cur.preds = pv.preds[e, 2:(length(pv.pred.species)+1)]

  
  
  #color bar
  library(s2dverification)
  
  svg(file = paste("cortex\\brain_size\\plasma_color_bar.svg",sep=""), width=8, height=1)
  
  ColorBar(brks=40, cols=color.scale, vertical=F, bar_limits=c(-0.75,0.75), subsampleg=6)
  
  dev.off()
  
  







dat = data.frame(Enhancer = enh.names[1:index-1], Pvalue = p.vals[1:index-1], Coeff = coeffs[1:index-1])
write.csv(dat, "liver\\pgls_lq_results.csv", row.names = FALSE)
write.csv(dat, sub(".csv", paste("_r", args[11], "_s", args[14], ".csv", sep=""),  args[10]), row.names = FALSE)





#Temoprary - report # of orthologs per species
#counts = double(length(pred.species))
#for (e in 1:nrow(preds)) {
#  cur.preds = preds[e, 2:(length(pred.species)+1)]
#  good.preds.indices = which(cur.preds != -1)
#  counts[good.preds.indices] = counts[good.preds.indices] + 1
#}
#dat = data.frame(Species=pred.species, Ortholog_Count=counts)
#write.csv(dat, "liver\\liver_ortholog_counts.csv", row.names=F)


