suppressMessages(library(ArchR))
ss <- function(x, pattern, slot = 1, ...) { sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("hg38")

###############################
### make an ArchR Projects ####
pd = read.delim('/projects/pfenninggroup/singleCell/BICCN_human_M1_SNARE-Seq2/rdas/Zhang_BICCN-H_20190523-20190611_huMOp_Final_Sample_Metadata.txt')
# subset to cells from final analyses
ArrowFiles = list.files(path = '/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/HumanDNase/M1SnareSeq/',pattern = '.arrow', full.names = TRUE)
  
# subset to cells from final analyses
proj = ArchRProject( ArrowFiles = ArrowFiles, outputDirectory = "ArchR_BICCN_huMOp_SNARE-Seq2", copyArrows = TRUE, showLogo = FALSE)
proj = saveArchRProject(ArchRProj = proj)
# 67,487 cells from replicate 1; 132,162 cells from replicate 2

proj = filterDoublets( ArchRProj = proj)
  
# subset to cells from final analyses
indKeep = match(pd$sample_name, toupper(ss(getCellNames(proj), '#',2)) )
indKeep = indKeep[!is.na(indKeep)]
  
# filter to cells annotated in final analyses 
proj = subsetArchRProject( ArchRProj = proj, cells = getCellNames(proj)[indKeep], outputDirectory = "ArchR_BICCN_huMOp_SNARE-Seq2_filtered")

# add cell metadata to ArchR project
pd = pd[match(toupper(ss(getCellNames(proj), '#',2)), pd$sample_name),]
for (name in names(pd)[-1]){
	proj = addCellColData(proj, data = pd[,name], name = name, cells = proj$cellNames, force = TRUE)
}
proj$subclass = make.names(proj$subclass)

proj = saveArchRProject(ArchRProj = proj)

# add iterative LSI
proj <- addIterativeLSI( ArchRProj = proj, useMatrix = "TileMatrix", 
                         name = "IterativeLSI",
                         LSIMethod = 2,
                         iterations = 4, # increase this if noticing subtle batch effects
                         scaleTo = 3000,
                         selectionMethod = 'var',
                         clusterParams = list( # See Seurat::FindClusters
                           resolution = c(0.2), # lower this if noticing subtle batch effects
                           sampleCells = 10000,  n.start = 10), 
                         varFeatures = 200000, # increase this if some sub-types are not separating
                         dimsToUse = 1:30, force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

# add harmony batch correction #
proj <- addHarmony( ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "patient", force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

# add umap
proj <- addUMAP( ArchRProj = proj, reducedDims = "Harmony", 
                 name = "UMAP", nNeighbors = 30, minDist = 0.5, 
                 metric = "cosine", force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

# add imputation
#proj <- addImputeWeights(proj, reducedDims = "Harmony")
#proj = saveArchRProject(ArchRProj = proj)

# make group coverage, call peaks, and 
# by subclass
proj<-addGroupCoverages(proj, groupBy="subclass", force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

# call peaks 
pdf('tmp.pdf')
proj<-addReproduciblePeakSet(proj, groupBy = "subclass", plot = FALSE)
dev.off()
proj = saveArchRProject(ArchRProj = proj)

# call peaks without grouping
proj = saveArchRProject(ArchRProj = proj, outputDirectory = "ArchR_BICCN_huMOp_SNARE-Seq2_filtered_noGroup")
proj<-addGroupCoverages(proj, force = TRUE)
pdf('tmp.pdf')
proj<-addReproduciblePeakSet(proj, plot = FALSE)
dev.off()
proj = saveArchRProject(ArchRProj = proj)

# add peak counts matrix 
#proj <- addPeakMatrix(proj)
#proj = saveArchRProject(ArchRProj = proj)

# add motif enrichment matrix
#proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
#proj = saveArchRProject(ArchRProj = proj)

# add motif deviations matrix
#proj <- addBgdPeaks(proj)
#proj <- addDeviationsMatrix(proj,  peakAnnotation = "Motif", force = TRUE)
#proj = saveArchRProject(ArchRProj = proj)

# add co-accessibility matrix
#proj <- addCoAccessibility(proj, reducedDims = "Harmony", dimsToUse = 1:30,
#                           scaleDims = TRUE, corCutOff = 0.75, k = 100, 
#                           knnIteration = 500, overlapCutoff = 0.8, 
#                           maxDist = 1e+05, scaleTo = 10^4, log2Norm = TRUE)
#proj = saveArchRProject(ArchRProj = proj)
