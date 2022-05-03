suppressMessages(library(ArchR))
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("mm10")

###############################
### make an ArchR Projects ####
ArrowFiles = list.files(path = '/projects/pfenninggroup/singleCell/BICCN_mouse_CATlas_snATAC-seq/arrow/', pattern = "2C|3C|4B|5D*.arrow", full.names = TRUE)

# subset to cells from final analyses
proj = ArchRProject( ArrowFiles = ArrowFiles, outputDirectory = "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/M1ScAtacSeqNEMO/ArchR_BICCN_mouse_NEMO_snATAC-Seq", copyArrows = TRUE, showLogo = FALSE)
proj = filterDoublets( ArchRProj = proj,  cutEnrich = .5)

# read in the metadata and filter
meta_table = read.table("/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/M1ScAtacSeqNEMO/CEMBA_MOp.L2.cluster.meta.txt", header=TRUE)
saveRDS(meta_table, "/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase/M1ScAtacSeqNEMO/CEMBA_MOp.L2.cluster.meta.RDS")
meta_fn =file.path('/projects/pfenninggroup/machineLearningForComputationalBiology/regElEvoGrant/MouseDNase', 'M1ScAtacSeqNEMO', 'CEMBA_MOp.L2.cluster.meta.RDS')
meta_df = readRDS(meta_fn)
meta_df$barcode = paste0(meta_df$Sample, "#", meta_df$barcode) # Do not need to filter for MOp because using only arrow files for MOp
meta_df$cluster = gsub('L23.a','L23', meta_df$cluster)
meta_df$cluster = gsub('L23.b','L23', meta_df$cluster)
meta_df$cluster = gsub('L23.c','L23', meta_df$cluster)
meta_df$cluster = gsub('IT.a','IT', meta_df$cluster)
meta_df$cluster = gsub('IT.b','IT', meta_df$cluster)
meta_df$cluster = gsub('_Arhgdib','', meta_df$cluster)
meta_df$cluster = gsub('_Mettl21e','', meta_df$cluster)
meta_df$cluster = gsub('_Ndnf','', meta_df$cluster)
meta_df$cluster = gsub('_Smad3','', meta_df$cluster)
meta_df$cluster = gsub('_Ntf3_Trim63','', meta_df$cluster)
meta_df$cluster = gsub('_Tac1','', meta_df$cluster)
meta_df$cluster = gsub('_Vsig2','', meta_df$cluster)
meta_df$cluster = gsub('_Chrna2_Myh8','', meta_df$cluster)
meta_df$cluster = gsub('_Man1a','', meta_df$cluster)
meta_df$cluster = gsub('_Stk33','', meta_df$cluster)
meta_df$cluster = gsub('_Chat','', meta_df$cluster)
meta_df$cluster = gsub('_Gcnt4','', meta_df$cluster)
meta_df$cluster = gsub('_Hcls1','', meta_df$cluster)
meta_df$cluster = gsub('_Lipg','', meta_df$cluster)

# rearrange cells
cellsKeep = intersect(proj$cellNames, meta_df$barcode)
proj = proj[cellsKeep,]
meta_df = meta_df[match(cellsKeep, meta_df$barcode),]
all.equal(meta_df$barcode, proj$cellNames)

for(i in seq(ncol(meta_df))){
  proj = addCellColData(
    ArchRProj = proj,
    data = meta_df[,i],
    cells = cellsKeep,
    name = names(meta_df)[i],
    force = T
  )
}

# add iterative LSI
# Requires 16GB
set.seed(1234)
proj <- addIterativeLSI(
  ArchRProj = proj, useMatrix = "TileMatrix",
  name = "IterativeLSI",
  LSIMethod = 2,
  iterations = 6, # increase this if noticing subtle batch effects
  scaleTo = 3000,
  selectionMethod = 'var',
  clusterParams = list( # See Seurat::FindClusters
    resolution = c(.1, .2, rep(.4, 3)), # lower this if noticing subtle batch effects
    sampleCells = 10000,  n.start = 10),
  varFeatures = 150000, # also can reduce this if noticing subtle batch effects
  dimsToUse = 1:40, force = TRUE)
  
# add harmony batch correction #
proj <- addHarmony( ArchRProj = proj, reducedDims = "IterativeLSI", name = "Harmony", groupBy = "Sample",force = TRUE)
					
# add imputation
#proj <- addImputeWeights(proj, reducedDims = "Harmony")

# add umap (swtich to ArchR 1.0.1 before doing this)
proj <- addUMAP( ArchRProj = proj, reducedDims = "Harmony", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)

# make group coverage
# Requires 8GB
proj<-addGroupCoverages(proj, groupBy="cluster", force = TRUE)
proj = saveArchRProject(ArchRProj = proj)

# add peak counts matrix
pdf('tmp.pdf')
proj<-addReproduciblePeakSet(proj, groupBy = "cluster", plot = FALSE)
dev.off()
proj = saveArchRProject(ArchRProj = proj)