##########################################
##### R program to create a database #####
##########################################


####### Versions ####################

# 1.1.1 - Basic version to compare anchors
# 1.1.2 - Extra function to test for vocal learning in brain chromatin
# 2.1.1 - Final bulk motor cortex predictions
# 2.1.3 - Add morgan's VL annotations
# 3.2.3 - Update for manuscript resubmission, t-test of differences, heatmap, coordinates
#enhPredVlCompare.2024.3.2.5clean.R - Adjust visualizations for final version, clean code

##########################################
##### Functions #####
##########################################

########### Gets breaks to use for heatmap #############
getBreaks <- function(min,max,scale){
	numVals <- (max - min)*(1/scale)
	x <- c(0:(numVals));
	y <- (x * scale) + min;
	return(y);
}


###### Makes a matrix ########
makeMat <- function(x,y) {
	mat <- matrix(NA,length(x),length(y));
	rownames(mat) <- x;
	colnames(mat) <- y;

	return(mat);
}

#########################################################
##### Parse Bed Input #####
#########################################################
bedFromRowFc <- function(rowNameV) {
	unlist(lapply(strsplit(gsub("_",".",rowNameV),"\\."),function(x) paste(x[2],":",x[3],"-",x[4],sep="")));
}



##########################################
##### Load appropriate packages #####
##########################################

library(gplots);
library(ggplot2);
library(GenomicRanges);
library(RColorBrewer);
library(VennDiagram);

library("ape");
library("phylolm");
library(scales)


#############################################################################################
##### Set directories #####
#############################################################################################

outFd <- "out32_vlRev/"
#workFd <- "work21/";

int2classV <- c("integer","numeric","character");

fileSuffix = "3.2.1";



#############################################################################################
##### Load the datasets #####
#############################################################################################

############### Load the neural cell types ##################

peakPred2L <- vector("list",5);
useCtV <- c("PV","mcx","liv");
names(peakPred2L) <- useCtV;


peakPredFn <- "/projects/pfenninggroup/mouseCxStr/NeuronSubtypeATAC/Zoonomia_CNN/predictions/PV_MoHu_AvgActivityPredictions_boreoeutheria.txt";
#peakPred2M <- as.matrix(read.delim(peakPredFn,row.names=1,header=T,stringsAsFactors=F,na.strings = "NA"));
peakPred2L[["PV"]] <-  as.matrix(read.delim(peakPredFn,row.names=1,header=T,stringsAsFactors=F,na.strings = "NA"));


lapply(peakPred2L,dim);

#Set the peak matrix to the default name
treeNames3V <- colnames(peakPred2L[["PV"]]);

#Get the common species across cell types
peakPred3L <- vector("list",0)
for(curCt in useCtV) {
	peakPred3L[[curCt]] <- peakPred2L[[curCt]][,treeNames3V];
}

############### Load the tissues ##################

tissuesV <- c("mcx");

#Level 1 = species of the coordinates
#Level 2 = species of in which the peaks originate
peakSpecV <- c("Homo_sapiens","Mus_musculus","Macaca_mulatta","Rousettus_aegyptiacus","Canis_lupus_familiaris");
peakInfoL <- vector("list",0);
peakInfoL[["mcx"]] <- vector("list",length(peakSpecV));
peakInfoL[["liv"]] <- vector("list",length(peakSpecV));
names(peakInfoL[["mcx"]]) <- peakSpecV;
names(peakInfoL[["liv"]] ) <- peakSpecV;

for(curTiss in tissuesV) {
	for(curSpec in peakSpecV) {
		peakInfoL[[curTiss]][[curSpec]] <- vector("list",0);
	}
}

refSpecL <- vector("list",0);
refSpecL[["mcx"]] <- c("mo","ma","ra","ba"); #Set of reference species to go through

refFnVL <- vector("list",0); #Where the reference file is with the species names
refFnVL[["mcx"]] <- c();
refFnVL[["mcx"]]["mo"] <- "m1Data/liftedPeakAllFileNamesConverted.txt";
refFnVL[["mcx"]]["ma"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt"
refFnVL[["mcx"]]["ra"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt"
refFnVL[["mcx"]]["ba"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/LiftedPeaksFixed/liftedPeakAllFileNamesConverted.txt"

replaceFnL <- vector("list",0);
#replaceFnL[["mcx"]] <- c("mo","ma","ra","ba"); #Set of reference species to go through
replaceFnL[["mcx"]]["mo"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/LiftedPeaks/Pfenning_bulk_Ctx_nonCDS_enhancerShort_";
replaceFnL[["mcx"]]["ma"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MacaqueAtac/OfM/peak/idr_reproducibility/LiftedPeaksFixed/idr.optimal_peak.inM1_nonCDS_enhancerShort_"
replaceFnL[["mcx"]]["ra"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/RatAtac/M1_AllReps/peak/idr_reproducibility/LiftedPeaksFixed/idr.optimal_peak_nonCDS_enhancerShort_"
replaceFnL[["mcx"]]["ba"] <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/BatAtac/OfMBat1K/call-reproducibility_idr/execution/LiftedPeaksFixed/idr.optimal_peak.inM1_nonCDS_enhancerShort_"


treeNames2L <- vector("list",0);
treeNames2L[["mcx"]] <-  vector("list",4); #Names of the trees for each reference species
names(treeNames2L[["mcx"]]) <- refSpecL[["mcx"]];

############ Read in mouse trees - Motor Cortex ###########################

#curRefSpec <- "ba";

for(curRefSpec in refSpecL[["mcx"]]) {

	#Read in the mouse tree names, fix
	treeNamesV <- as.vector(as.matrix(read.delim(refFnVL[["mcx"]][curRefSpec],stringsAsFactors=F,header=F)));
	treeNames2V <- gsub(replaceFnL[["mcx"]][curRefSpec],"",treeNamesV);
	treeNames2V <- gsub("_summitExtendedMin50Max2XProtect5_GenBankNames.bed.gz","",treeNames2V);
	treeNames2V <- gsub("_summitExtendedMin50Max2XProtect5.bed.gz","",treeNames2V);
	treeNames2V <- gsub("_summitExtendedMin50Max2XProtect5_RefSeqNames.bed.gz","",treeNames2V);
	treeNames2V <- gsub("_summitExtendedMin50Max2XProtect5_UCSCNames.bed.gz","",treeNames2V);

	#Fix macaque tree names
	treeNames2V <- gsub("GenBankNames_","",treeNames2V);

	#Fix bat tree names
	treeNames2V <- gsub("200MBat_SequenceNames_","",treeNames2V);


	treeNames2V[1] <- "Mus_musculus";
	treeNames2V[119] <- "Macaca_mulatta";
	treeNames2V[201] <- "Rattus_norvegicus";
	treeNames2V[205] <- "Rousettus_aegyptiacus";


	treeNames2L[["mcx"]][[curRefSpec]] <- treeNames2V;
	#treeNames2V

	for(curPeakSpec in peakSpecV) {
		curPeakInfoFn <- treeNamesV[which(treeNames2V == curPeakSpec)];
		curPeakInfoF <- read.delim(curPeakInfoFn,stringsAsFactors=F,header=F);
		colnames(curPeakInfoF) <- c("chr","start","stop","peakId","na1","na2","na3","na4","na5","bedScore");
		curPeakInfoF$peakId <- paste(curPeakInfoF$peakId,curRefSpec,sep="_");
		rownames(curPeakInfoF) <- curPeakInfoF$peakId;

		peakInfoL[["mcx"]][[curPeakSpec]][[curRefSpec]] <- curPeakInfoF;

	}
}

#peakInfoL[["Homo_sapiens"]][["mo"]][1:5,]
#peakInfoL[["Mus_musculus"]][["mo"]][1:5,]
lapply(peakInfoL[["mcx"]],function(x) lapply(x,head));

#Collapse peak info into one data frame per species
peakInfo2L <- vector("list",0);
peakInfo2L[["mcx"]] <- vector("list",0);

for(curPeakSpec in peakSpecV) {
	peakInfo2L[["mcx"]][[curPeakSpec]] <- data.frame();
	for(curRefSpec in names(peakInfoL[["mcx"]][[curPeakSpec]])) {
		peakInfo2L[["mcx"]][[curPeakSpec]] <- rbind(peakInfo2L[["mcx"]][[curPeakSpec]],peakInfoL[["mcx"]][[curPeakSpec]][[curRefSpec]]);
	}
}


#Read in the matrix Information
peakPred1L <- vector("list",0);
peakPred1L[["mcx"]] <- as.matrix(read.delim("m1Data/MoMaRaBaEnhancersNRSummit.txt",header=F,stringsAsFactors=F,row.names=1,na.strings=c("NA","-1")));

tmpPeakFn <- "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/Cortex_All_ATAC_out/peak/macs2/idr/optimal_set/Pfenning_bulk_Ctx_nonCDS_enhancerShort_peakPredictionsMultiSpecies.txt";
tmpM1PeakPredM <- as.matrix(read.delim(tmpPeakFn,header=F,stringsAsFactors=F,row.names=1,na.strings=c("NA","-1")))

colnames(peakPred1L[["mcx"]]) <-  treeNames2L[["mcx"]][["mo"]][match(peakPred1L[["mcx"]]["peak16505_mo",],tmpM1PeakPredM["peak16505",])];


colnames(peakPred1L[["mcx"]])[which(colnames(peakPred1L[["mcx"]]) == "Canis_lupus_familiaris")] <- "Canis_lupus";

peakPred3L[["mcx"]] <- peakPred1L[["mcx"]][,treeNames3V]


########### Genome qual stats ###########
genomeQualStatsFn <- "/home/apfennin/projects/covidTest/orthGene/zoonomia/basicGenomeQualStats.2.1.1.v2simp.csv"
genomeQualStatsF <- read.csv(genomeQualStatsFn,stringsAsFactors=F,header=T)
genomeQualStats2F <- genomeQualStatsF[match(treeNames3V,genomeQualStatsF$Species),]; #Indexed by treeNamesV

genomeQualStats3F <- genomeQualStats2F;
rownames(genomeQualStats3F) <- genomeQualStats3F$Species

#grep("canis",genomeQualStats2F$Species)

############  Species evolutionary distance ############
specEvolDistFn <- "/home/apfennin/projects/covidTest/orthGene/zoonomia/200Mammal_hal_tree_noancestors_matrix.csv"
specEvolDistM <- as.matrix(read.csv(specEvolDistFn,stringsAsFactors=F,header=T,row.names=1));
specEvolDist2M <- specEvolDistM;
matrixSpecV <- rownames(specEvolDistM);
for(curSpec1 in 1:(length(matrixSpecV)-1)) {
	for(curSpec2 in (curSpec1+1):length(matrixSpecV)) {
		specEvolDist2M[curSpec1,curSpec2] <- specEvolDist2M[curSpec2,curSpec1]
	}
}


######## Load species information ############
######## Load Morgan's vocal learning annotations ############

vlAnnotF <- read.csv("data/Boreoeutheria_VLtrait_MW.csv")

genomeTraitF <- genomeQualStats2F[,c("Species","Common.Name","Clade")];
genomeTraitF$vl <- vlAnnotF[match(genomeTraitF$Species,gsub(" ","_",vlAnnotF$Species.Name)),"Vocal_Learner"]

genomeTrait2F <- genomeTraitF;
rownames(genomeTrait2F) <- genomeTrait2F$Species

#sort(table(genomeTrait2F$Species))

######## Load Clade-specific vocal learing information, from Morgan's code ############

hillConvertFn <- "/projects/pfenninggroup/vocalLearning/4exp_RER/data/Zoonomia-Hiller_species_name_conversions.csv";
hillConvertF <- read.csv(hillConvertFn,stringsAsFactors=F);
hillConvertF$SpeciesId <- gsub(" ","_",hillConvertF$Species.binomials);

if(F) {
	match(treeNames3V,hillConvertF$SpeciesId);
	setdiff(treeNames3V,hillConvertF$SpeciesId);

	genomeTraitF[match("Platanista_gangetica",genomeTraitF$Species),]

}

#Manually Convert Mismatching names

hillConvertF[which(hillConvertF$SpeciesId == "Bison_bison_bison"),"SpeciesId"] <- "Bison_bison";
hillConvertF[which(hillConvertF$SpeciesId == "Balaenoptera_acutorostrata_scammoni"),"SpeciesId"] <- "Balaenoptera_acutorostrata";
hillConvertF[which(hillConvertF$SpeciesId == "Canis_lupus_familiaris"),"SpeciesId"] <- "Canis_lupus";
hillConvertF[which(hillConvertF$SpeciesId == "Cebus_capucinus_imitator"),"SpeciesId"] <- "Cebus_capucinus";
hillConvertF[which(hillConvertF$SpeciesId == "Ceratotherium_simum_simum"),"SpeciesId"] <- "Ceratotherium_simum";
hillConvertF[which(hillConvertF$SpeciesId == "Colobus_angolensis_palliatus"),"SpeciesId"] <- "Colobus_angolensis";
hillConvertF[which(hillConvertF$SpeciesId == "Dicerorhinus_sumatrensis_sumatrensis"),"SpeciesId"] <- "Dicerorhinus_sumatrensis";
hillConvertF[which(hillConvertF$SpeciesId == "Enhydra_lutris_nereis"),"SpeciesId"] <- "Enhydra_lutris";
hillConvertF[which(hillConvertF$SpeciesId == "Gorilla_gorilla_gorilla"),"SpeciesId"] <- "Gorilla_gorilla";
hillConvertF[which(hillConvertF$SpeciesId == "Marmota_marmota_marmota"),"SpeciesId"] <- "Marmota_marmota";
hillConvertF[which(hillConvertF$SpeciesId == "Murina_aurata_feae"),"SpeciesId"] <- "Murina_feae";
hillConvertF[which(hillConvertF$SpeciesId == "Neophocaena_asiaeorientalis_asiaeorientalis"),"SpeciesId"] <- "Neophocaena_asiaeorientalis";
hillConvertF[which(hillConvertF$SpeciesId == "Panthera_tigris_altaica"),"SpeciesId"] <- "Panthera_tigris";
hillConvertF[which(hillConvertF$SpeciesId == "Perognathus_longimembris_pacificus"),"SpeciesId"] <- "Perognathus_longimembris";
hillConvertF[which(hillConvertF$SpeciesId == "Peromyscus_maniculatus_bairdii"),"SpeciesId"] <- "Peromyscus_maniculatus";
hillConvertF[which(hillConvertF$SpeciesId == "Peromyscus_maniculatus_bairdii"),"SpeciesId"] <- "Peromyscus_maniculatus";

genomeTrait2F$vlAnalysisStatus <- genomeTrait2F$Clade;

genomeTrait2F[which(genomeTraitF$Species == "Homo_sapiens"),"vlAnalysisStatus"] <- "vl_human";

#genomeQualStats3F[names(sort(specEvolDist2M["Rousettus_aegyptiacus",]))[1:60],c("Species","Common.Name")]

batsV <- c("Rousettus_aegyptiacus","Pteropus_alecto","Eidolon_helvum","Pteropus_vampyrus","Macroglossus_sobrinus","Hipposideros_armiger",
"Craseonycteris_thonglongyai","Rhinolophus_sinicus","Hipposideros_galeritus","Megaderma_lyra","Noctilio_leporinus",
"Micronycteris_hirsuta","Anoura_caudifer","Carollia_perspicillata","Eptesicus_fuscus","Lasiurus_borealis",
"Desmodus_rotundus","Tonatia_saurophila","Myotis_davidii","Artibeus_jamaicensis","Tadarida_brasiliensis",
"Pteronotus_parnellii","Miniopterus_schreibersii","Miniopterus_natalensis","Myotis_myotis","Myotis_lucifugus",
"Murina_feae","Mormoops_blainvillei","Pipistrellus_pipistrellus","Myotis_brandtii");

genomeTrait2F[batsV,"vlAnalysisStatus"] <- "vl_bats";

#genomeQualStats3F[names(sort(specEvolDist2M["Zalophus_californianus",]))[1:30],c("Species","Common.Name")]
sealsV <- c("Odobenus_rosmarus","Neomonachus_schauinslandi","Mirounga_angustirostris","Leptonychotes_weddellii","Zalophus_californianus");
genomeTrait2F[sealsV,"vlAnalysisStatus"] <- "vl_seals";

#genomeQualStats3F[names(sort(specEvolDist2M["Eschrichtius_robustus",]))[1:30],c("Species","Common.Name")]
whalesV <- c("Balaenoptera_bonaerensis","Balaenoptera_acutorostrata","Eubalaena_japonica","Kogia_breviceps",
"Orcinus_orca","Monodon_monoceros","Delphinapterus_leucas","Neophocaena_asiaeorientalis","Phocoena_phocoena",
"Inia_geoffrensis","Lipotes_vexillifer","Ziphius_cavirostris","Mesoplodon_bidens","Platanista_gangetica",
"Tursiops_truncatus");
genomeTrait2F[whalesV,"vlAnalysisStatus"] <- "vl_whales";

table(genomeTrait2F$vlAnalysisStatus);

genomeTrait2F[which(is.na(genomeTrait2F$vl)),"vlAnalysisStatus"] <- "exclude";

table(genomeTrait2F$vlAnalysisStatus);

######## Load Clade-specific vocal learing information, from Morgan's code ############

#treeUrl <- "http://www.andrew.cmu.edu/user/apfennin/2021_Fall_GEB/hw2/Zoonomia_ChrX_lessGC40_241species_30Consensus.tree";
#system(paste("wget ",treeUrl));
zoonomiaTree <- read.tree(file = "Zoonomia_ChrX_lessGC40_241species_30Consensus.tree") #read tree in Newick format
zoonomiaTree

if(F) {
	#loadedDataWsFn <- paste("work/","loadedData.3.2.1.1.ws",sep=""); #The basic normalized data
	#loadedDataWsFn <- paste("work/","loadedData.3.2.1.2.ws",sep=""); #The basic normalized data, fixed bug with misannoted VL species
	loadedDataWsFn <- paste("work/","loadedData.3.2.3.1.ws",sep=""); #The basic normalized data, fixed names bug
	#save.image(loadedDataWsFn);
	load(loadedDataWsFn);


}


#############################################################################################
##### Create basic functions to use for the comparison #####
#############################################################################################

############## Functions used to compute differences ################

#Compuates a t-test by a single individual to another group
#testThresh - must as have at least that many species to compute
#Computes refSpec - species spec
tTestFromMat <- function(specSetOutV,refSpec,testThresh) {
	specDiffV <- refSpec - specSetOutV;
	specDiffV <- specDiffV[which(!is.na(specDiffV))];

	if(length(specDiffV) >= testThresh && var(specDiffV) > 0) {
		return(t.test(specDiffV));
	} else
		return(NA);

}

#test is mean x - mean y
#Computes specSetTestV - specSetOutV which is reference - outgroup
#Compuates a t-test by a test group individual to another group
#testThresh - must as have at least that many species to compute

tTestFromMatMult <- function(specSetTestV,specSetOutV,testThreshRef,testThreshOut) {
	specSetTestSubV <- specSetTestV[which(!is.na(specSetTestV))];
	specSetOutSubV <- specSetOutV[which(!is.na(specSetOutV))];

	if(length(specSetTestSubV) >= testThreshRef && length(specSetOutSubV) >= testThreshOut) {
		return(t.test(x=specSetTestSubV,y=specSetOutSubV));
	} else
		return(NA);

}

#autmatically creates a dataframe of species differences
#specNameOutV - a list of species names for the outgrup
#refSpecName - The reference species name
createDataDiffSingle <- function(specNameOutV,refSpecName,testThresh,fnPredM) {
	#Compute the p-value differences
	curTestL <- apply(fnPredM,1,function(x) tTestFromMat(x[specNameOutV],x[refSpecName],testThresh));

	#Extract the p-values
	curTestV <- unlist(lapply(curTestL,function(x) if(!is.na(x)) {x$p.value} else {NA}))
	curTest2V <- curTestV[which(!is.na(curTestV))];

	#Create a dataframe
	curResF <- data.frame(peakId=names(curTest2V),pval=curTest2V,padj=p.adjust(curTest2V,method="fdr"),stringsAsFactors=F);
	curResF$specPred <- fnPredM[curResF$peakId,refSpecName];
	curResF$outPred <- apply(fnPredM[curResF$peakId,specNameOutV],1,mean,na.rm=T);
	curResF$grpDiff <- curResF$specPred - curResF$outPred;

	#Add coordinates
	if(F) {
		curCoordM <- makeMat(curResF$peakId,peakSpecV);
		for(curPeakSpec in peakSpecV) {
			tmpPeakF <- peakInfo2L[[curPeakSpec]]
			curCoordM[,curPeakSpec] <-  paste(tmpPeakF[curResF$peakId,"chr"],":",tmpPeakF[curResF$peakId,"start"],"-",tmpPeakF[curResF$peakId,"stop"],sep="");
		}
		colnames(curCoordM) <- paste("coord",colnames(curCoordM),sep="_");

		curResF <- cbind(curResF,curCoordM);
	}

	return(curResF);

}

#autmatically creates a dataframe of species differences
#specNameOutV - a list of species names for the outgrup
#refSpecName - The reference species name
createDataDiffMult<- function(specNameRefV,specNameOutV,testThreshRef,testThreshOut,fnPredM) {
	#Compute the p-value differences
	curTestL <- apply(fnPredM,1,function(x) tTestFromMatMult(x[specNameRefV], x[specNameOutV],testThreshRef,testThreshOut));

	#Extract the p-values
	curTestV <- unlist(lapply(curTestL,function(x) if(!is.na(x)) {x$p.value} else {NA}))
	curTest2V <- curTestV[which(!is.na(curTestV))];

	#Create a dataframe
	curResF <- data.frame(peakId=names(curTest2V),pval=curTest2V,padj=p.adjust(curTest2V,method="fdr"),stringsAsFactors=F);
	curResF$specPred <- apply(fnPredM[curResF$peakId,specNameRefV],1,mean,na.rm=T);
	curResF$outPred <- apply(fnPredM[curResF$peakId,specNameOutV],1,mean,na.rm=T);
	curResF$grpDiff <- curResF$specPred - curResF$outPred;

	#Add coordinates
	if(F) {
		curCoordM <- makeMat(curResF$peakId,peakSpecV);
		for(curPeakSpec in peakSpecV) {
			tmpPeakF <- peakInfo2L[[curPeakSpec]]
			curCoordM[,curPeakSpec] <-  paste(tmpPeakF[curResF$peakId,"chr"],":",tmpPeakF[curResF$peakId,"start"],"-",tmpPeakF[curResF$peakId,"stop"],sep="");
		}
		colnames(curCoordM) <- paste("coord",colnames(curCoordM),sep="_");

		curResF <- cbind(curResF,curCoordM);
	}

	return(curResF);

}

######### Process the trait info ###############

genomeTrait3F <- genomeTrait2F;

genomeTrait3F$curVlTrait <- genomeTrait3F$vl;

vlAnalysisLevelsV <- c("Euarchonta","Glire","Laurasiatheria","vl_bats","vl_whales","vl_seals","vl_human","exclude");
genomeTrait3F$vlAnalysisStatusT <- factor(genomeTrait3F$vlAnalysisStatus,levels=vlAnalysisLevelsV);


######################### Identify differences across each clade ########################
############## Annotate motor cortex for specific clades ################

finalTissuesV <- c("mcx","PV");

diffTtestFrameL <- vector("list",0); #List of data frames with full test results

for(curPredName in finalTissuesV) {

	diffTtestFrameL[[curPredName]] <- vector("list",0);

	#curPredName <- "mcx";
	#curPredM <- peakPred3L[[curPredName]][,colnames(peakPred3L[["mcx"]])]
	curPredM <- peakPred3L[[curPredName]];

	#table(genomeTrait3F$vlAnalysisStatus)

	#sort(specEvolDist2M["Homo_sapiens",])

	### Annotate human vs. Euarchonta ###
	diffTtestFrameL[[curPredName]][["vl_human"]] <- createDataDiffSingle(genomeTrait3F[which(genomeTrait3F$vlAnalysisStatus == "Euarchonta"),"Species"],"Homo_sapiens",3,curPredM);
	#quantile(diffResFrameL[["vl_human"]]$grpDiff,c(0:10)*.1)

	laurasV <- genomeTrait3F[which(genomeTrait3F$vlAnalysisStatus == "Laurasiatheria"),"Species"]
	diffTtestFrameL[[curPredName]][["vl_bats"]] <- createDataDiffMult(genomeTrait3F[which(genomeTrait3F$vlAnalysisStatus == "vl_bats"),"Species"],laurasV,3,3,curPredM);
	diffTtestFrameL[[curPredName]][["vl_whales"]] <- createDataDiffMult(genomeTrait3F[which(genomeTrait3F$vlAnalysisStatus == "vl_whales"),"Species"],laurasV,3,3,curPredM);
	diffTtestFrameL[[curPredName]][["vl_seals"]] <- createDataDiffMult(genomeTrait3F[which(genomeTrait3F$vlAnalysisStatus == "vl_seals"),"Species"],laurasV,3,3,curPredM);

	#lapply(diffResFrameL,function(x) length(which(x$padj < 0.01 & x$grpDiff > 0)));
	#lapply(diffResFrameL,function(x) length(which(x$padj < 0.01 & x$grpDiff < 0)));

}


#Print tables for publication
if(F) {
	for(curPredName in finalTissuesV) {
		for(curCompName in names(diffTtestFrameL[[curPredName]])){
			curDiffFn <- paste(outFd,"table/tacitTtest.",curPredName,".",curCompName,".1.csv",sep="");
			write.csv(diffTtestFrameL[[curPredName]][[curCompName]],file=curDiffFn);
		}
	}

}


#############################################################################################
##### Study VL phylolm output #####
#############################################################################################
#Daniel and Irene's results, with published way of adjusting p-value


################# Load in the data ################

m1PeakPermFn <- "data/cortex_vl_results_combined_mousefix_12-31_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner_adjp.csv";
pvPeakPermFn <- "data/pv_vl_results_combined_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner_adjp.csv"

permulationsResL <- vector("list",0);

permulationsResL[["mcx"]] <- read.csv(m1PeakPermFn,stringsAsFactors=F);
permulationsResL[["PV"]] <- read.csv(pvPeakPermFn,stringsAsFactors=F);

ttestPvalThresh <- 0.05;
ttestQvalThresh <- 1;

permPvalThresh <- 1;
permQvalThresh <- 0.1;

################# Figure S7 - Plot value distributions ################

curTissue <- "mcx";
curPvalHistF <- permulationsResL[[curTissue]];

curPvalHistF$associationDirection <- "pos";
curPvalHistF[which(curPvalHistF$Coeff < 0),"associationDirection"] <- "neg"


if(F) {


	pdf(paste(outFd,"pvalHist/hist.postperm.phyloGlmP.",curTissue,".",fileSuffix,".3vlfix.pdf",sep=""),width=5,height=4);
		curP <- ggplot(curPvalHistF, aes(x = Exp_Pvalue, fill=associationDirection)) +
		#geom_histogram(bins=50, color=curPvalHistF$associationDirection, position="dodge") +
		geom_histogram(color="black",position='identity',alpha=0.5) +
		#geom_density(adjust=0.5,fill="black") +
		theme_bw();
		#scale_x_continuous(breaks=c(25,50,75,100),limits=c(20, 105));
		print(curP);
	dev.off();

}



################# Create Venn Diagramns about which regions are in which species ################

#m1PeakPerm2F <- m1PeakPermF[order(m1PeakPermF$bh),]
#curPredName <- "mcx"
#curCategory <- "vl_bats";

vennResultsL <- vector("list",0);

for(curPredName in names(permulationsResL)) {

	#curCategoryV <- names(diffTtestFrameL[[curPredName]]);
	curCategoryV <- c("vl_bats","vl_whales","vl_seals","vl_human")

	vennResultsL[[curPredName]] <- vector("list",0);

	#Create lists of things that are higher in vocal learners in each category
	curUpPeakV <- permulationsResL[[curPredName]][which(permulationsResL[[curPredName]]$Coeff > 0 & permulationsResL[[curPredName]]$bh <= permQvalThresh),"Enhancer"];
	curUpPeakL <- vector("list",0); #For each vocal learning clade, does the peak reach significance
	for(curCategory in curCategoryV) {
		tmpTtestResF <- diffTtestFrameL[[curPredName]][[curCategory]][curUpPeakV,]
		curUpPeakL[[curCategory]] <- tmpTtestResF[which(tmpTtestResF$grpDiff > 0 & tmpTtestResF$pval <= ttestPvalThresh),"peakId"]
	}

	curDownPeakV <- permulationsResL[[curPredName]][which(permulationsResL[[curPredName]]$Coeff < 0 & permulationsResL[[curPredName]]$bh <= permQvalThresh),"Enhancer"];
	curDownPeakL <- vector("list",0); #For each vocal learning clade, does the peak reach significance
	for(curCategory in curCategoryV) {
		tmpTtestResF <- diffTtestFrameL[[curPredName]][[curCategory]][curDownPeakV,]
		curDownPeakL[[curCategory]] <- tmpTtestResF[which(tmpTtestResF$grpDiff < 0 & tmpTtestResF$pval <= ttestPvalThresh),"peakId"]
	}

		vennResultsL[[curPredName]][["up"]] <- curUpPeakL;
		vennResultsL[[curPredName]][["down"]] <- curDownPeakL;

}

directionsV <- c("up","down");

for(curPredName in names(permulationsResL)) {

	for(curDirection in directionsV) {

		if(F) {

			currVennL <- vennResultsL[[curPredName]][[curDirection]];
			currVennFn <- paste(outFd,"venn/vennDiagram.",curPredName,".",curDirection,".",fileSuffix,".2reord.png",sep="");
			venn.diagram(x=currVennL,
					category.names = names(currVennL),
					fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
					filename = currVennFn,
					imagetype ="png"

			);
		}
	#	dev.off();

	}

}

#### Merge Venn results for mcx and PV ####

vennResultsAllUpL <- vector("list",0);
vennResultsAllDownL <- vector("list",0);

for(curClade in names(vennResultsL$PV$up)) {
	vennResultsAllUpL[[curClade]] <- c(vennResultsL[["mcx"]][["up"]][[curClade]],vennResultsL[["PV"]][["up"]][[curClade]]);
	vennResultsAllDownL[[curClade]] <- c(vennResultsL[["mcx"]][["down"]][[curClade]],vennResultsL[["PV"]][["down"]][[curClade]]);

}

vennResultsMergeL <- vector("list",0);
vennResultsMergeL[["up"]] <- vennResultsAllUpL;
vennResultsMergeL[["down"]] <-vennResultsAllDownL;

for(curDirection in names(vennResultsMergeL)) {

	if(T) {

		currVennL <- vennResultsMergeL[[curDirection]];
		currVennFn <- paste(outFd,"venn/vennDiagram.merge.",curDirection,".",fileSuffix,".2reord.png",sep="");
		venn.diagram(x=currVennL,
				height = 1500, width = 2000,
				category.names = names(currVennL),
				fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
				filename = currVennFn,
				imagetype ="png"

		);
	#	dev.off();

	}

}

################# Print out Mouse Bed Files of Peaks for Annotate ################

mouseBedL <- vector("list",0);
mouseBedL[["mcx"]] <- peakInfo2L[["mcx"]][["Mus_musculus"]][permulationsResL[["mcx"]]$Enhancer,]
mouseBedL[["mcx"]] <- mouseBedL[["mcx"]][which(!is.na(mouseBedL[["mcx"]][,1])),];
#lapply(peakInfo2L[["mcx"]],dim)
mouseBedL[["mcx"]][1:5,]

curPeakFn <- paste("intermedFiles/postProcPeak.","mcx",".","mm10",".1.bed",sep="")

if(F) {
	write.table(mouseBedL[["mcx"]],file=curPeakFn,row.names=F,col.names=F,sep="\t",quote=F);
}

tmpPvEnh <- gsub("\\.","_",permulationsResL[["PV"]]$Enhancer)
mouseBedL[["PV"]] <- data.frame(do.call(rbind,strsplit(tmpPvEnh,split=c("\\_")))[,2:4])
colnames(mouseBedL[["PV"]]) <- c("chr","start","stop")
mouseBedL[["PV"]]$peakId <- permulationsResL[["PV"]]$Enhancer
mouseBedL[["PV"]] <- mouseBedL[["PV"]][grep("mm10",mouseBedL[["PV"]]$peakId),]

curPeakFn <- paste("intermedFiles/postProcPeak.","pv",".","mm10",".1.bed",sep="")

if(F) {
	write.table(mouseBedL[["PV"]],file=curPeakFn,row.names=F,col.names=F,sep="\t",quote=F);
}


################# Print out Mouse Bed Files of Peaks for Annotate ################

#export PATH=$PATH:/home/apfennin/bin/homer/bin/

##awk -F '\t' '{print $2}' * | sort | uniq -c | sort -nr
#awk -vOFS='\t' -vFS='\t' '{print $1,$2,$3,$4}' intermedFiles/postProcPeak.mcx.mm10.1.bed > intermedFiles/postProcPeak.mcx.mm10.2simp.bed


#/home/apfennin/bin/homer/bin/annotatePeaks.pl intermedFiles/postProcPeak.mcx.mm10.2simp.bed mm10 -annStats intermedFiles/geneAnt/mcx.permPeakAnt.1.stats.txt -gtf /home/apfennin//tools/mm10_ant/Mus_musculus.GRCm38.87.2.gtf  > intermedFiles/geneAnt/mcx.permPeakAnt.1.txt
#/home/apfennin/bin/homer/bin/annotatePeaks.pl intermedFiles/postProcPeak.mcx.mm10.2simp.bed mm10 -annStats intermedFiles/geneAnt/mcx.permPeakAnt.2prot.stats.txt -gtf /home/apfennin//tools/mm10_ant/Mus_musculus.GRCm38.87.2.prot.gtf  > intermedFiles/geneAnt/mcx.permPeakAnt.2prot.txt
#/home/apfennin/bin/homer/bin/annotatePeaks.pl intermedFiles/postProcPeak.pv.mm10.1.bed mm10 -annStats intermedFiles/geneAnt/pv.permPeakAnt.2prot.stats.txt -gtf /home/apfennin//tools/mm10_ant/Mus_musculus.GRCm38.87.2.prot.gtf  > intermedFiles/geneAnt/pv.permPeakAnt.2prot.txt


######### Gene anntoations - Protein Only - MCX ##############
mcxAntFn <- "intermedFiles/geneAnt/mcx.permPeakAnt.2prot.txt";
#mcxAntF <- read.delim(mcxAntFn,stringsAsFactors=F,header=F,na.strings=c(".","NA"),skip=1,colClasses=int2classV[c(3,3,3,1,3,3,1,3,3,3,1,3,3,1)]);
mcxAntF <- read.delim(mcxAntFn,stringsAsFactors=F,header=T,na.strings=c(".","NA"));
colnames(mcxAntF)[1] <- "peakId";
#colnames(mcxAntF) <- c("mid", "tssEns", "tssSym", "tssDist","tssProtEns", "tssProtSym", "tssProtDist", "loc", "fullEns", "fullSym", "fullDist", "fullProtEns","fullProtSym","fullProtDist");
#motAntF$mid <- gsub("/home/japostol/projects/motoshi/peaks/","",motAntF$mid)
#rownames(mcxAntF) <- motAntF$mid

transcript2symF <- read.csv("data/ensembl.transcript2sym.GRCm39.csv",stringsAsFactors=F)

mcxAnt2F <- cbind(mcxAntF,transcript2symF[match(mcxAntF$Nearest.PromoterID, transcript2symF$Transcript.stable.ID),]);
rownames(mcxAnt2F) <- mcxAnt2F$peakId

#Rename peaks based on nearest mouse gene
uGenesV <- unique(mcxAnt2F$Gene.name);

numGeneV <- rep(0,length(uGenesV));
names(numGeneV) <- uGenesV;

mcxAnt2F$peakGeneId <- mcxAnt2F$peakId;

for(curPeakNum in 1:length(mcxAnt2F$peakId)) {
	curGene <- mcxAnt2F[curPeakNum,"Gene.name"];
	curGeneNum <- numGeneV[curGene] + 1;
	numGeneV[curGene] <- curGeneNum;
	mcxAnt2F[curPeakNum,"peakGeneId"] <- paste(curGene,curGeneNum,sep="_");
}

mcxAnt2F[1:5,]

mcxResAntF <-  cbind(permulationsResL[["mcx"]][match(mcxAnt2F$peakId,permulationsResL[["mcx"]]$Enhancer),],mcxAnt2F);

mcxResAnt2permordF <- mcxResAntF[order(mcxResAntF$Exp_Pvalue),]; #order based on permulations p

mcxResAnt3singleGeneF <- mcxResAnt2permordF[match(uGenesV,mcxResAnt2permordF$Gene.name),] ;#only keep the most signficant enhancer per gene
mcxResAnt3singleGeneF$permScore <- (1 - mcxResAnt3singleGeneF$Exp_Pvalue) * sign(mcxResAnt3singleGeneF$Coeff);

mcxResAnt4reordF <- mcxResAnt3singleGeneF[order(mcxResAnt3singleGeneF$permScore),]

mcxResAnt5reordF <- mcxResAnt3singleGeneF[order(mcxResAnt3singleGeneF$Exp_Pvalue,decreasing=F),];

mcxResAnt6reordTmpF <- mcxResAnt5reordF[which(mcxResAnt5reordF$bh <= 0.1),];

mcxResAnt6reordTmpF[,c("Gene.Name","bh")]


if(F) {
	curPeakFn <- paste("intermedFiles/postProcPeak.","mcx",".","mm10",".topSig.",1000,".1.bed",sep="")
	write.table(mouseBedL[["mcx"]][match(mcxResAnt5reordF[1:1000,"Enhancer"],mouseBedL[["mcx"]]$peakId),c(1:4)],file=curPeakFn,row.names=F,col.names=F,sep="\t",quote=F);

}


######### Gene anntoations - Protein Only - PV ##############
pvAntFn <- "intermedFiles/geneAnt/pv.permPeakAnt.2prot.txt";
#mcxAntF <- read.delim(mcxAntFn,stringsAsFactors=F,header=F,na.strings=c(".","NA"),skip=1,colClasses=int2classV[c(3,3,3,1,3,3,1,3,3,3,1,3,3,1)]);
pvAntF <- read.delim(pvAntFn,stringsAsFactors=F,header=T,na.strings=c(".","NA"));
colnames(pvAntF)[1] <- "peakId";
#colnames(mcxAntF) <- c("mid", "tssEns", "tssSym", "tssDist","tssProtEns", "tssProtSym", "tssProtDist", "loc", "fullEns", "fullSym", "fullDist", "fullProtEns","fullProtSym","fullProtDist");
#motAntF$mid <- gsub("/home/japostol/projects/motoshi/peaks/","",motAntF$mid)
#rownames(mcxAntF) <- motAntF$mid

transcript2symF <- read.csv("data/ensembl.transcript2sym.GRCm39.csv",stringsAsFactors=F)

pvAnt2F <- cbind(pvAntF,transcript2symF[match(pvAntF$Nearest.PromoterID, transcript2symF$Transcript.stable.ID),]);
rownames(pvAnt2F) <- pvAnt2F$peakId

#Rename peaks based on nearest mouse gene
uGenesV <- unique(pvAnt2F$Gene.name);

numGeneV <- rep(0,length(uGenesV));
names(numGeneV) <- uGenesV;

pvAnt2F$peakGeneId <- pvAnt2F$peakId;

for(curPeakNum in 1:length(pvAnt2F$peakId)) {
	curGene <- pvAnt2F[curPeakNum,"Gene.name"];
	curGeneNum <- numGeneV[curGene] + 1;
	numGeneV[curGene] <- curGeneNum;
	pvAnt2F[curPeakNum,"peakGeneId"] <- paste(curGene,curGeneNum,sep="_");
}

pvAnt2F[1:5,]

pvAnt2F <-  cbind(permulationsResL[["PV"]][match(pvAnt2F$peakId,permulationsResL[["PV"]]$Enhancer),],pvAnt2F);

pvResAnt2permordF <- pvAnt2F[order(pvAnt2F$Exp_Pvalue),]; #order based on permulations p


############### Print out table of regulatory elements #######################

fullMcxResF <- permulationsResL[["mcx"]][which(permulationsResL[["mcx"]]$Pvalue < 1 & permulationsResL[["mcx"]]$bh < 0.1),]
fullMcxResF$Trials <- NA;
fullMcxResF <- cbind(fullMcxResF,mcxAnt2F[match(fullMcxResF$Enhancer,mcxAnt2F$peakId),c("Chr","Start","End","Distance.to.TSS","Gene.name","Gene.stable.ID")]);

if(F) {
	curResFn <- paste(outFd,"table/tacit.sigPermResults.","mcx",".","mm10.2.csv",sep="");
	write.csv(fullMcxResF,file=curResFn,row.names=T,col.names=T,quote=F);

}

fullPvResF <- permulationsResL[["PV"]][which(permulationsResL[["PV"]]$Pvalue < 1 & permulationsResL[["PV"]]$bh < 0.1),]
fullPvResF$Trials <- NA;
fullPvResF <- cbind(fullPvResF,pvAnt2F[match(fullPvResF$Enhancer,pvAnt2F$peakId),c("Chr","Start","End","Distance.to.TSS","Gene.name","Gene.stable.ID")]);

if(F) {
	curResFn <- paste(outFd,"table/tacit.sigPermResults.","pv",".","mm10.2.csv",sep="")
	write.csv(fullPvResF,file=curResFn,row.names=T,col.names=T,quote=F);

}


#############################################################################################
##### Create Heatmap of Relevant Regulatory Elements #####
#############################################################################################

curTissue <- "mcx";

#Get the predictions for the right peaks, species order
curRawPredM <- peakPred3L[[curTissue]];

#Confirm that nothing is in the tree that isn't in the matrix;
removeSpeciesV <- setdiff(zoonomiaTree$tip.label,colnames(curRawPredM));

#Subset the tree with the species for which there are predictions
zoonomiaSubTree <- drop.tip(zoonomiaTree,removeSpeciesV);

curSpeciesOrderV <- zoonomiaSubTree$tip.label;


#Get the set of peaks to use
curDownPeakV <- permulationsResL[[curTissue]][which(permulationsResL[[curTissue]]$Coeff < 0 & permulationsResL[[curTissue]]$bh <= permQvalThresh),"Enhancer"];
curUpPeakV <- permulationsResL[[curTissue]][which(permulationsResL[[curTissue]]$Coeff > 0 & permulationsResL[[curTissue]]$bh <= permQvalThresh),"Enhancer"];

curPlotPeaksV <- c(curDownPeakV,curUpPeakV);

curRawPred2subM <- curRawPredM[curPlotPeaksV,curSpeciesOrderV]; #Re-order based on peak height and phylogenetic tree

curGenomeTraitF <- genomeTrait2F[curSpeciesOrderV,]

table(genomeTrait2F$vl)

vlColorsV <- rep(NA,length(curSpeciesOrderV));
names(vlColorsV) <- curSpeciesOrderV;
vlColorsV[which(curGenomeTraitF$vl == 1)] <- "red";
vlColorsV[which(curGenomeTraitF$vl == 0)] <- "black";
vlColorsV[which(is.na(curGenomeTraitF$vl))] <- "gray";

curEnhColorsV <- c(rep("purple",length(curDownPeakV)),rep("green",length(curUpPeakV)));

#Normalize relative to vocal non-learners
curNormPred2subM <- curRawPred2subM - apply(curRawPred2subM[,which(curGenomeTraitF$vl == 0)],1,mean,na.rm=T)

#Rename rows to correspond to nearest mouse gene
curNormPred3renameM <- curNormPred2subM;
newRowNamesV <- mcxResAnt3singleGeneF[match(rownames(curNormPred2subM),mcxResAnt3singleGeneF$Enhancer),"peakGeneId"];
rownames(curNormPred3renameM)[which(!is.na(newRowNamesV))] <- newRowNamesV[which(!is.na(newRowNamesV))]
rownames(curNormPred3renameM)[37] <- rownames(curNormPred2subM)[37]

curNormPred4scaleM <- curNormPred3renameM*(1/apply(curNormPred3renameM,1,function(x) max(abs(x),na.rm=T)))


#### Print out heatmap  ####
if(F) {

	#Heatmap with mean normalized for vocal non-learners
	pdf(paste(outFd,"differentialRegionHeatmap.meanNorm.mcx.",fileSuffix,".2prot.pdf",sep=""),width=6,height=6);

		heatmap.2(x=curNormPred3renameM , Rowv=T, Colv=F, scale="none", dendrogram="none",
			symm = F,
			RowSideColors=curEnhColorsV,
			ColSideColors=vlColorsV,
			#col="bluered",breaks=getBreaks(0,1,0.01), na.color="gray",
			col="bluered",breaks=getBreaks(-1,1,0.01), na.color="gray",
			distfun = function(dat) as.dist(1-cor(t(dat),use='pairwise.complete.obs')),

			trace="none",density.info="none",
			);
	dev.off();

	#apply(curNormPred3renameM,1,function(x) max(abs(x),na.rm=T))

}

#############################################################################################
##### Enrichment of motor cortex open chromatin #####
#############################################################################################

curTissue <- "mcx"

######## Get the set of peaks to use - Motor Cortex ############

cellSectPThresh <- 0.05
cellsectDownPeakV <- permulationsResL[[curTissue]][which(permulationsResL[[curTissue]]$Coeff < 0 & permulationsResL[[curTissue]]$Exp_Pvalue <= cellSectPThresh),"Enhancer"];
cellsectUpPeakV <- permulationsResL[[curTissue]][which(permulationsResL[[curTissue]]$Coeff > 0 & permulationsResL[[curTissue]]$Exp_Pvalue <= cellSectPThresh),"Enhancer"];


cellsectDownPeakV <- intersect(cellsectDownPeakV,mouseBedL[["mcx"]]$peakId)
cellsectUpPeakV <- intersect(cellsectUpPeakV,mouseBedL[["mcx"]]$peakId)

cellsectNonPeakV <- setdiff(intersect(permulationsResL[[curTissue]]$Enhancer,mouseBedL[["mcx"]]$peakId),c(cellsectDownPeakV,cellsectUpPeakV));

mcxPeakSetsL <- vector("list",0);
mcxPeakSetsL[["vl_down"]] <- cellsectDownPeakV;
mcxPeakSetsL[["vl_up"]] <- cellsectUpPeakV;
mcxPeakSetsL[["vl_non"]] <- cellsectNonPeakV;
mcxPeakSetsL[["vl_down_rand"]] <- sample(cellsectNonPeakV,length(cellsectDownPeakV),replace=F);

mxcPeakNamesV <- names(mcxPeakSetsL);

mouseMcxPeakGrL <- vector("list",0);

for(curMcx in names(mcxPeakSetsL)) {
	tmpPeakBedF <- mouseBedL[["mcx"]][match(mcxPeakSetsL[[curMcx]],mouseBedL[["mcx"]]$peakId),]
	mouseMcxPeakGrL[[curMcx]] <-  with(tmpPeakBedF, GRanges(chr, IRanges(start+1, stop), id=peakId))
}
#lapply(mouseMcxPeakGrL,head);

######## Get the set of peaks to use - Single Cell Motor Cortex ############

mouseSnPeakFnV <- system("ls /projects/pfenninggroup/machineLearningForComputationalBiology/Cortex_Cell-TACIT/data/raw_data/filtered_peaks_final/BICCN_Mouse/*",intern=T)
mouseSnCellTypeV <- gsub("/projects/pfenninggroup/machineLearningForComputationalBiology/Cortex_Cell-TACIT/data/raw_data/filtered_peaks_final/BICCN_Mouse/BICCN_mouse_Cortex_snATAC.","",mouseSnPeakFnV)

mouseSnCellTypeV <- gsub(".mm10.enhancer.narrowPeak.gz","",mouseSnCellTypeV)


names(mouseSnPeakFnV) <- mouseSnCellTypeV;


mouseCxSnrnaSeqPeakGrL <- vector("list",0);

for(curCellType in mouseSnCellTypeV) {
	curFn <- mouseSnPeakFnV[curCellType];
	tmpCellTypeNarrowPeakF <- read.delim(curFn,header=F,stringsAsFactors=F);
	mouseCxSnrnaSeqPeakGrL[[curCellType]] <-  with(tmpCellTypeNarrowPeakF, GRanges(V1, IRanges(V2+1, V3), id=V4))

}

lapply(mouseCxSnrnaSeqPeakGrL,length);

####### Intersect the motor cortex peaks with the mouse single cell peaks ############

peakSectNumM <- makeMat(mxcPeakNamesV,mouseSnCellTypeV);
peakPropNumM <- peakSectNumM;

numPerm <- 1000;

for(curMcx in mxcPeakNamesV) {
	for(curCellType in mouseSnCellTypeV) {
		tmpSectV <- findOverlaps(	mouseMcxPeakGrL[[curMcx]],mouseCxSnrnaSeqPeakGrL[[curCellType]],ignore.strand=T,select="first",type="any");
		#tmpSectV <- findOverlaps(mouseCxSnrnaSeqPeakGrL[[curCellType]],mouseMcxPeakGrL[[curMcx]],ignore.strand=T,select="first");
		tmpSectV <- tmpSectV[which(!is.na(tmpSectV))];
		tmpSectV <- unique(tmpSectV);

		#tmpSectV <- countOverlaps(	mouseMcxPeakGrL[[curMcx]],mouseCxSnrnaSeqPeakGrL[[curCellType]],ignore.strand=T,type="any");
		#peakSectNumM[curMcx,curCellType] <- length(which(tmpSectV > 0));
		peakSectNumM[curMcx,curCellType] <- length(tmpSectV);
		peakPropNumM[curMcx,curCellType] <- peakSectNumM[curMcx,curCellType]/length(mouseCxSnrnaSeqPeakGrL[[curCellType]]);
	}
}
peakSectNumM

#Generate Null Distribution for decreasing peaks
tempRandDecM <- makeMat(1:numPerm,mouseSnCellTypeV);
for(curPerm in 1:numPerm) {
	tmpRandPeakV <- sample(cellsectNonPeakV,length(cellsectDownPeakV),replace=F);
	tmpRandPeakF <- mouseBedL[["mcx"]][match(tmpRandPeakV,mouseBedL[["mcx"]]$peakId),]
	tmpRandPeakGr <-  with(tmpRandPeakF, GRanges(chr, IRanges(start+1, stop), id=peakId))

	for(curCellType in mouseSnCellTypeV) {
		tmpSectV <- findOverlaps(	tmpRandPeakGr,mouseCxSnrnaSeqPeakGrL[[curCellType]],ignore.strand=T,select="first");
		#tmpSectV <- findOverlaps(mouseCxSnrnaSeqPeakGrL[[curCellType]],mouseMcxPeakGrL[[curMcx]],ignore.strand=T,select="first");
		tmpSectV <- tmpSectV[which(!is.na(tmpSectV))];
		tmpSectV <- unique(tmpSectV);
		tempRandDecM[curPerm,curCellType] <- length(tmpSectV);
		#peakPropNumM[curMcx,curCellType] <- peakSectNumM[curMcx,curCellType]/length(mouseCxSnrnaSeqPeakGrL[[curCellType]]);
	}
}

#Generate Null Distribution for increasing peaks
tempRandIncM <- makeMat(1:numPerm,mouseSnCellTypeV);
for(curPerm in 1:numPerm) {
	tmpRandPeakV <- sample(cellsectNonPeakV,length(cellsectUpPeakV),replace=F);
	tmpRandPeakF <- mouseBedL[["mcx"]][match(tmpRandPeakV,mouseBedL[["mcx"]]$peakId),]
	tmpRandPeakGr <-  with(tmpRandPeakF, GRanges(chr, IRanges(start+1, stop), id=peakId))

	for(curCellType in mouseSnCellTypeV) {
		tmpSectV <- findOverlaps(	tmpRandPeakGr,mouseCxSnrnaSeqPeakGrL[[curCellType]],ignore.strand=T,select="first");
		#tmpSectV <- findOverlaps(mouseCxSnrnaSeqPeakGrL[[curCellType]],mouseMcxPeakGrL[[curMcx]],ignore.strand=T,select="first");
		tmpSectV <- tmpSectV[which(!is.na(tmpSectV))];
		tmpSectV <- unique(tmpSectV);
		tempRandIncM[curPerm,curCellType] <- length(tmpSectV);
		#peakPropNumM[curMcx,curCellType] <- peakSectNumM[curMcx,curCellType]/length(mouseCxSnrnaSeqPeakGrL[[curCellType]]);
	}
}

decPermPV <- c();
incPermPV <- c();
for(curCellType in mouseSnCellTypeV) {
	decPermPV[curCellType] <- length(which(peakSectNumM["vl_down",curCellType] <= tempRandDecM[,curCellType]))
	incPermPV[curCellType] <- length(which(peakSectNumM["vl_up",curCellType] <= tempRandIncM[,curCellType]))
}

sort(decPermPV)/numPerm;
sort(incPermPV)/numPerm;


if(F) {


	#wsPermuteResFn <- paste("work/","permuate.allVl.3.2.3.1.ws",sep=""); #Processes data load
	wsPermuteResFn <- paste("work/","permuate.allVl.3.2.3.2p05.ws",sep=""); #Processes data load
	#save.image(wsPermuteResFn);
	load(wsPermuteResFn);


}

# Old - original barplots#
cellTypeResF <- data.frame(cellType = mouseSnCellTypeV,
														vl_down_num = peakSectNumM["vl_down",],
														vl_up_num = peakSectNumM["vl_up",]
													)

cellTypeResF$enrich_down <- log2(peakSectNumM["vl_down",]/apply(tempRandDecM,2,mean,na.rm=T));
cellTypeResF$enrich_up <- log2(peakSectNumM["vl_up",]/apply(tempRandIncM,2,mean,na.rm=T));
cellTypeResF$permP_down <- decPermPV/1000;
cellTypeResF$permP_up <- incPermPV/1000;
cellTypeResF$cellTypeT <- factor(cellTypeResF$cellType,levels=names(sort(decPermPV)));
cellTypeResF$sigDownP <- "notSig";
cellTypeResF$sigUpP <- "notSig";
cellTypeResF[which(cellTypeResF$permP_down < 0.05),"sigDownP"] <- "p<0.05";
cellTypeResF[which(cellTypeResF$permP_down < 0.01),"sigDownP"] <- "p<0.01";
cellTypeResF[which(cellTypeResF$permP_up < 0.05),"sigUpP"] <- "p<0.05";
cellTypeResF[which(cellTypeResF$permP_up < 0.01),"sigUpP"] <- "p<0.01";

cellTypeResF$sigDownPT <- factor(cellTypeResF$sigDownP,levels=c("notSig","p<0.05","p<0.01"))
cellTypeResF$sigUpPT <- factor(cellTypeResF$sigUpP,levels=c("notSig","p<0.05","p<0.01"))

#write.csv(cellTypeResF,paste(outFd,"table/cellTypeEnrichments.",fileSuffix,".1.csv",sep=""));

if(F) {

	pdf(paste(outFd,"cellType/cellTypeEnrich.vl_decreasing",fileSuffix,".3p05.pdf",sep=""),width=5,height=3);
		curP <- ggplot(cellTypeResF, aes(x = cellTypeT , y = enrich_down, fill=sigDownPT)) +
		geom_bar(stat='identity') +
		#geom_smooth(method = "lm", se=FALSE, color="blue") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90,hjust=1)) +
		scale_fill_manual(values=c("gray","#997570","red"));
		#scale_fill_gradient2(low="red",	mid = "grey50",high = "blue",	midpoint = 0.5,	na.value = "grey50", guide = "colourbar");
		#geom_line(data=curKnlData2F, aes(x=age, y=pred), colour="red", size=1.5);
		print(curP);
	dev.off();

	pdf(paste(outFd,"cellType/cellTypeEnrich.vl_increasing",fileSuffix,".3p05.pdf",sep=""),width=5,height=3);
		curP <- ggplot(cellTypeResF, aes(x = cellTypeT , y = enrich_up, fill=sigUpPT)) +
		geom_bar(stat='identity') +
		#geom_smooth(method = "lm", se=FALSE, color="blue") +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90,hjust=1)) +
		scale_fill_manual(values=c("gray","#997570","red"));
		#geom_line(data=curKnlData2F, aes(x=age, y=pred), colour="red", size=1.5);
		print(curP);
	dev.off();


}

#### Barplots with error bars requested for 2024 review - Add replicate information ####

fullPermDataIncF <- data.frame();
fullPermDataDecF <- data.frame();

for(curCellType in mouseSnCellTypeV) {

	tmpPermIncF <- data.frame(cellType=rep(curCellType,numPerm),permHits=tempRandIncM[,curCellType]);
	tmpPermIncF$enrich <- cellTypeResF[curCellType,"enrich_up"];
	tmpPermIncF$perm <- cellTypeResF[curCellType,"perm_up"];
	tmpPermIncF$sigP <- cellTypeResF[curCellType,"sigUpP"];
	tmpPermIncF$realHits <- cellTypeResF[curCellType,"vl_up_num"];
	fullPermDataIncF <- rbind(fullPermDataIncF,tmpPermIncF);

	tmpPermDecF <- data.frame(cellType=rep(curCellType,numPerm),permHits=tempRandDecM[,curCellType]);
	tmpPermDecF$enrich <- cellTypeResF[curCellType,"enrich_down"];
	tmpPermDecF$perm <- cellTypeResF[curCellType,"perm_down"];
	tmpPermDecF$sigP <- cellTypeResF[curCellType,"sigDownP"];
	tmpPermDecF$realHits <- cellTypeResF[curCellType,"vl_down_num"];
	fullPermDataDecF <- rbind(fullPermDataDecF,tmpPermDecF)

}

fullPermDataIncF$cellTypeT <- factor(fullPermDataIncF$cellType,levels=names(sort(decPermPV,decreasing=T)));
fullPermDataDecF$cellTypeT <- factor(fullPermDataDecF$cellType,levels=names(sort(decPermPV,decreasing=T)));

if(F) {

	pdf(paste(outFd,"cellType/cellTypeEnrich.errorBars2024.vl_decreasing.",fileSuffix,".3p05.1.pdf",sep=""),width=8,height=3);
		curP <- ggplot(fullPermDataDecF, aes(x = cellTypeT , y = permHits)) +
		geom_violin() +
		#geom_bar(stat='identity') +
		#geom_smooth(method = "lm", se=FALSE, color="blue") +
		#stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
		geom_boxplot(width=0.1) +
		geom_point(aes(x=cellTypeT, y=realHits),color='red') +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90,hjust=1));
		#scale_fill_manual(values=c("gray","#997570","red"));
		#scale_fill_gradient2(low="red",	mid = "grey50",high = "blue",	midpoint = 0.5,	na.value = "grey50", guide = "colourbar");
		#geom_line(data=curKnlData2F, aes(x=age, y=pred), colour="red", size=1.5);
		print(curP);
	dev.off();

	pdf(paste(outFd,"cellType/cellTypeEnrich.errorBars2024.vl_increasing.",fileSuffix,".3p05.1.pdf",sep=""),width=8,height=3);
		curP <- ggplot(fullPermDataIncF, aes(x = cellTypeT , y = permHits)) +
		geom_violin() +
		#geom_bar(stat='identity') +
		#geom_smooth(method = "lm", se=FALSE, color="blue") +
		#stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
		geom_boxplot(width=0.1) +
		geom_point(aes(x=cellTypeT, y=realHits),color='red') +
		theme_bw() +
		theme(axis.text.x = element_text(angle = 90,hjust=1));
		#scale_fill_manual(values=c("gray","#997570","red"));
		#scale_fill_gradient2(low="red",	mid = "grey50",high = "blue",	midpoint = 0.5,	na.value = "grey50", guide = "colourbar");
		#geom_line(data=curKnlData2F, aes(x=age, y=pred), colour="red", size=1.5);
		print(curP);
	dev.off();

}
