suppressMessages(library(ArchR))
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("hg38")

proj = loadArchRProject("ArchR_BICCN_huMOp_SNARE-Seq2_filtered")

#########################
# extract IDR peak set
cellTypes = levels(factor(proj$subclass))
cellTypesCrossSpeciesNeuron = c("L2.3.IT", "L5.ET", "L5.IT", "L5.6.NP", "L6b", "L6.CT", "L6.IT", "L6.IT.Car3", "LAMP5", "PVALB", "SNCG", "SST", "SST.CHODL", "VIP")
idrPeakFiles = file.path('ArchR_BICCN_huMOp_SNARE-Seq2_filtered','PeakCalls',paste0(cellTypesCrossSpeciesNeuron,'-reproduciblePeaks.gr.rds'))
names(idrPeakFiles) = cellTypesCrossSpeciesNeuron
idrPeakList = GRangesList(lapply(idrPeakFiles, readRDS))

# write IDR peaks to bed file
for (name in names(idrPeakList)){
	bedFileName = file.path(paste0('BICCN_huMOp_SNARE-Seq2_',name,'_reproduciblePeaks.bed'))
	rtracklayer::export(idrPeakList[[name]], bedFileName)
}
