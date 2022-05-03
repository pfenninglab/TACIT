suppressMessages(library(ArchR))
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("mm10")

proj = loadArchRProject("ArchR_BICCN_mouse_NEMO_snATAC-Seq")

#########################
# extract IDR peak set
cellTypes = levels(factor(proj$cluster))
cellTypesCrossSpeciesNeuron = c("Chodl", "L23", "L5.IT", "L5.PT", "L6.CT", "L6.IT", "Lamp5", "NP", "Pv", "Sncg", "Sst", "Vip")
idrPeakFiles = file.path('ArchR_BICCN_mouse_NEMO_snATAC-Seq','PeakCalls',paste0(cellTypesCrossSpeciesNeuron,'-reproduciblePeaks.gr.rds'))
names(idrPeakFiles) = cellTypesCrossSpeciesNeuron
idrPeakList = GRangesList(lapply(idrPeakFiles, readRDS))

# write IDR peaks to bed file
for (name in names(idrPeakList)){
	bedFileName = file.path(paste0('BICCN_mouse_NEMO_snATAC-Seq_',name,'_reproduciblePeaks.bed'))
	rtracklayer::export(idrPeakList[[name]], bedFileName)
}