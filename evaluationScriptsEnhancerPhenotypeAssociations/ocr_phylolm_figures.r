#Contains code used for producing figures showing OCR-phenotype associations
#Includes hard-coded file names, which need to be changed for each tissue and phenotype
#The example used is cortex and brain size, unless specified.

library(ape)
library(stringr)
library(phylolm)
library(viridis)
library(s2dverification)

#Read tree
tree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")


#Read phenotype - exaple: brain size
traits = read.csv("zoonomia_phenotypes_4-16-21.csv")
trait.col= "Brain.resid"
trait.all = traits[[trait.col]]
iv = which(trait.all != "NA") 
trait = trait.all[iv]
species.lower = traits$species.binomial[iv]
trait.species = str_to_title(species.lower)
names(trait) = trait.species
#Collect taxnomy data
clades.all = traits$Taxonomic.lineage
names(clades.all) = str_to_title(traits$species.binomial)


#Read OCR activity predictions
preds = read.csv(file = "cortex\\MoMaRaBaEnhancersNRSummit.txt", header = F, sep = "\t")
names(preds)[1] = "Enhancer"
te = read.csv(file = "BoreoeutheriaSpeciesNamesNew.txt", header=F)
pred.species = gsub(" ", "_", te$V1)
names(preds)[2:(length(pred.species)+1)] = pred.species

common.species = intersect(intersect(pred.species, tree$tip.label), trait.species)
te = which(trait.species %in% common.species)
tree.common = keep.tip(tree, common.species)



#SCATTER PLOT (example: cortex / brain size)
e = 75902 #Specific row of predictions containing enhancer of interest
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
cet = which(grepl("Cetacea", int.clades))
shape[ape] = 15
shape[cet] = 17

laur.line = phylolm(Y ~ X, data = dat[laur,], phy=int.tree.di, model = "BM")
if (length(glir) > 0) {
glir.line = phylolm(Y ~ X, data = dat[glir,], phy=int.tree.di, model = "BM")
}
euar.line = phylolm(Y ~ X, data = dat[euar,], phy=int.tree.di, model = "BM")

svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_v3.svg",sep=""), width=8, height=8)
#change ylim depending on trait - for brain size, -0.75 - 0.75
plot(dat$X, dat$Y, col=color, pch=shape, xlim=c(0, 1), ylim=c(-0.75,0.75), xlab="OCR activity prediction", cex=2, ylab = "Brain size residual", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
abline(a=m$coefficients[1], b=m$coefficients[2], lwd=4)
abline(a=laur.line$coefficients[1], b=laur.line$coefficients[2], lty=2, col="#ff8000", lwd=2)
if (length(glir) > 0) {
abline(a=glir.line$coefficients[1], b=glir.line$coefficients[2], lty=2, col="#0000ff", lwd=2)
}
abline(a=euar.line$coefficients[1], b=euar.line$coefficients[2], lty=2, col="#008000", lwd=2)
#legend("bottom", x.intersp=1, xjust=0, yjust=0, bty="n", text.width=c(0.14,0.05,0.12), legend=c("Laurasiatheria", "Glires", "Euarchonta"), col=c("red", "blue", "darkgreen"), pch=c(19,19,19), text.col="black", horiz=T)
#rect(0.2,-0.9,0.8,-0.65)

dev.off()
  
  

#PV SCATTER (also use this to load alternate-format PV data for figures below)
#Read PV predictions & phenotypes (this is for an alternate format)
pv.table = read.csv("PVenhPreds_Euarchontoglires_brainResidSig.csv")
pv.enh = "mm10_chr13.114793237.114793737" #Name of specific PV OCR
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

if (length(glir) > 0) {
  glir.line = phylolm(Y ~ X, data = dat[glir,], phy=int.tree.di, model = "BM")
}
euar.line = phylolm(Y ~ X, data = dat[euar,], phy=int.tree.di, model = "BM")

svg(file = paste("cortex\\brain_size\\pv_brain-size_",name,"_v3.svg",sep=""), width=8, height=8)
#change ylim depending on trait - for brain size, -0.75 - 0.75
plot(dat$X, dat$Y, col=color, pch=shape, xlim=c(0, 1), ylim=c(-0.75,0.75), xlab="OCR activity prediction", cex=2, ylab = "Brain size residual", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
abline(a=m$coefficients[1], b=m$coefficients[2], lwd=4)
if (length(glir) > 0) {
  abline(a=glir.line$coefficients[1], b=glir.line$coefficients[2], lty=2, col="#0000ff", lwd=2)
}
abline(a=euar.line$coefficients[1], b=euar.line$coefficients[2], lty=2, col="#008000", lwd=2)
dev.off()



#NEW-SYTLE FIGURE (prediction vs clade with phenotype coloring)
e = 75902 #Row corresponding to specific OCR (cortex/liver)
print(e)
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

dev.off()



#COLOR BAR

svg(file = paste("cortex\\brain_size\\plasma_color_bar.svg",sep=""), width=8, height=1)
ColorBar(brks=40, cols=color.scale, vertical=F, bar_limits=c(-0.75,0.75), subsampleg=6)
dev.off()



#INSET SCATTER PLOTS
whales = grep("Cetacea", int.clades)
if (length(dat$X[whales]) > 0){
svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_whales.svg",sep=""), width=6, height=6)
#To get actual inset, used ylim=c(0.25,0.8)
plot(dat$Y[whales], dat$X[whales], col=colors[whales], pch=19, ylim=c(0, 1), xlim=c(-0.75,0.75),xlab="Brain size residual", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
line = lm(X ~ Y, dat[whales,])
abline(a=line$coefficients[1], b=line$coefficients[2], lwd=4)
dev.off()
}

apes = grep("Hominoidea", int.clades)
svg(file = paste("cortex\\brain_size\\cortex_brain-size_",name,"_apes.svg",sep=""), width=6, height=6)
plot(dat$Y[apes], dat$X[apes], col=colors[apes], pch=19, ylim=c(0, 1), xlim=c(-0.75,0.75),xlab="Brain size residual", cex=3, ylab = "OCR Activity Prediction", main=name, cex.lab=1.5, cex.main = 2, cex.axis=2)
line = lm(X ~ Y, dat[apes,])
abline(a=line$coefficients[1], b=line$coefficients[2], lwd=4)
dev.off()
  


#SOLITARY PLOTS
#Read phenotype
soltraits = read.csv("cortex\\solitary\\SocialPhenotypes.csv")
iv = which(soltraits$solitary_living_euarchontoglires != "NA")
soltrait = soltraits$solitary_living_euarchontoglires[iv]
names(soltrait) = soltraits$Species.Name[iv]
common.species = intersect(intersect(pred.species, tree$tip.label), names(soltrait))
tree.common = keep.tip(tree, common.species)

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

int.clades = clades.all[rownames(dat)]


#Hardcode a few missing clades (cortex)
names(int.clades)[9] = "Capromys_pilorides"
int.clades[9] = "Hystricomorpha"
int.clades[24] = "Hystricomorpha"
names(int.clades)[24] = "Cuniculus_paca"
names(int.clades)[77] = "Tupaia_tana"
int.clades[77] = "Scandentia"

#(pv)
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

dev.off()


#Load PV stuff for above
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




