#Script for generating a tree of species given pairwise distances. 

library(ape)
library(stringr)


#Parse pairwise distances for all species
mat = read.csv("cortex_species_cosine_dist.csv", header=F)
d = as.dist(mat)
t = hclust(d, method="average")
sp = read.csv("cortex_species.txt", header=F, sep="\t")
t$labels = sp$V2

species = gsub(" ", "_", sp$V2)

cladeList = c("Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Tylopoda", "Suina", "Ruminantia", "Cetacea", "Megachiroptera", "Microchiroptera", "Lagomorpha", "Sciuromorpha", "Hystricomorpha", "Castorimorpha", "Myomorpha", "Strepsirrhini", "Platyrrhini", "Cercopithecoidea", "Hominoidea", "Scandentia")
cladeNums = c(1:20) #c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 12, 13, 14, 15, 16)
#g = 14
# = c(1, g, 2, 3, g, g, 4, 5, 6, g, 7, 8, g, 9, 10, 11, 12, 13)

#Parse clades from phenotype data file
traits = read.csv("..\\zoonomia_phenotypes_4-16-21.csv")
clades.all = traits$Taxonomic.lineage
names(clades.all) = str_to_title(traits$species.binomial)

int.clades = clades.all[species]

names(int.clades)[31] = "Capromys_pilorides"
int.clades[31] = "Hystricomorpha"
names(int.clades)[52] = "Crocidura_indochinensis"
int.clades[52] = "Eulipotyphla"
names(int.clades)[56] = "Cuniculus_paca"
int.clades[56] = "Hystricomorpha"
names(int.clades)[133] = "Murina_feae"
int.clades[133] = "Chiroptera;Microchiroptera"
names(int.clades)[178] = "Platanista_gangetica"
int.clades[178] = "Cetacea"
names(int.clades)[213] = "Tupaia_tana"
int.clades[213] = "Scandentia"

cladeX = rep(0, 222)
for (i in 1:222) {
  cladeNum = cladeNums[str_detect(int.clades[i], cladeList)] #uses stringr
  if (length(cladeNum) != 1) {
    print(int.clades[i])
  }
  cladeX[i] = cladeNum
}

#rb = rainbow(19)
#rb = c(rainbow(13), "darkgray")
rb = c("gray50", "goldenrod", "sienna", "deeppink", "darkorange", "tomato", "darkred", "red", "magenta", "darkviolet", "darkcyan", "cyan", "dodgerblue", "blue", "darkblue", "limegreen", "mediumseagreen", "green", "darkgreen", "gold4")


svg(file = paste("cortex_dendrogram_full-color.svg",sep=""), width=6, height=11)
#par(mar=c(5.1, 4.1, 4.1, 2.1), mgp=c(50, 40, 0))
#plot(as.dendrogram(t), horiz=T, axes=F, xlab="", ylab="", sub="", main = "Species dendrogram inferred from cortex OCR activity predictions")
plot(as.phylo(t), tip.color=rb[cladeX], cex = 0.4, no.margin=T)

dev.off()
