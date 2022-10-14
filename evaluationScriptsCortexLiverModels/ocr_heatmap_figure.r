library(ape)
library(stringr)

tree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree")


##For reading zoonomia sheet
traits = read.csv("zoonomia_phenotypes_4-16-21.csv")
clades.all = traits$Taxonomic.lineage
names(clades.all) = str_to_title(traits$species.binomial)

te = read.csv(file = "BoreoeutheriaTreeNamesNew.txt", header=F)
pred.species = gsub(" ", "_", te$V1)
preds = read.csv("cortex\\cortex_sample_1k_ward_preds.csv", sep=",", header=F)
names(preds) = pred.species


#             5,              2            28,           6                2           1        10              11          20           3             5               11              3                14           10               11             15                  5
cladeList = c("Eulipotyphla", "Pholidota", "Carnivora", "Perissodactyla", "Tylopoda", "Suina", "Ruminantia", "Cetacea", "Chiroptera", "Lagomorpha", "Sciuromorpha", "Hystricomorpha", "Castorimorpha", "Myomorpha", "Strepsirrhini", "Platyrrhini", "Cercopithecoidea", "Hominoidea", "Scandentia")

cladeNums = c(1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 10, 11, 11, 12, 13, 14, 15, 16)

clades.all["Capromys_pilorides"] = "Hystricomorpha"
clades.all["Crocidura_indochinensis"] = "Eulipotyphla"
clades.all["Cuniculus_paca"] = "Hystricomorpha"
clades.all["Murina_feae"] = "Chiroptera"    
clades.all["Tupaia_tana"] = "Scandentia"
clades.all["Platanista_gangetica"] = "Cetacea"

int.clades = clades.all[pred.species]
cladeX = rep(0, length(pred.species))
for (i in 1:length(pred.species)) {
  cladeNum = cladeNums[str_detect(int.clades[i], cladeList)] #uses stringr
  if (length(cladeNum) != 1) {
    print(i)
  } else {
    cladeX[i] = cladeNum
  }
}

#We want to preserve intra-clade order while obtaining new ordering of clades
#This is kinda gross but gets the job done
new.order = c()
new.cladeX = c()
for (i in 1:16) {
  new.order = c(new.order, pred.species[which(cladeX == i)])
  new.cladeX = c(new.cladeX, cladeX[which(cladeX == i)])
}
new.preds = preds[,new.order]

#Now, identify points where the clade changes
col.lab = rep("", length(new.cladeX))
for (i in 2:length(new.cladeX)) {
  if (new.cladeX[i] != new.cladeX[i-1]) {
    col.lab[i] = "*"
  }
}

temp1 <- as.matrix(new.preds)
#rand <- sample(nrow(temp1))
colors = c("#f2f2f2", "#ebdbda", "#e4c4c3", "#dcadad", "#d39696", "#c98081", "#be696c", "#b35258", "#a63944", "#991b32")
my_palette <-append(rep("#808080", 10), colors)

#Using gplots 3.0.1 from github/ChristophH/gplots, which fixes a bug that makes heatmap.2 really slow
library("gplots")

#small size
png(file = "cortex\\cortex_sample_1k_new5.png", width=2625, height=2040)
heatmap.2(temp1, col=my_palette, dendrogram="none", trace='none', symm=F,symkey=F,symbreaks=T, scale="none", key=F, labRow=F, labCol=col.lab, cexCol=1,
          Colv=F, Rowv=F, margins=c(10,4), lhei = c(1,15), lwid=c(1,40))
          
          #margins=c(4,4), lhei = c(1,40), lwid=c(1,40))#, margins=c(10,4), lhei = c(1,12), lwid=c(1,20))
#title("All species", line=-2.25, cex.main = 9)
#mtext("222 Boreoeutherians", line=3.5, side=1, cex=6)
#mtext("1000 cortex OCRs", line=-0.5, side=2, cex=6)
#rect(0.0178, 0.0134, 1.0012, 0.9596, lwd=40, ljoin=1, lend=2, xpd=T)
dev.off()
