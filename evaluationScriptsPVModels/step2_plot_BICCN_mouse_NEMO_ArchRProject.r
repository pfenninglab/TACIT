suppressMessages(library(ArchR))
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("mm10")

proj = loadArchRProject("ArchR_BICCN_mouse_NEMO_snATAC-Seq")

#########################################
# 1) heatmap of cell subclass by experiment_short
library(pheatmap)
cM <- confusionMatrix(paste0(proj$cluster), paste0(proj$Sample))

pdf('BICCN_mouse_NEMO_snATAC-Seq_combined_cluster_unnomralized_sample_heatmap.pdf', width = 7, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM)), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

cM <- cM / Matrix::rowSums(cM)

pdf('BICCN_mouse_NEMO_snATAC-Seq_combined_cluster_sample_heatmap.pdf', width = 7, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM)), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

pdf('BICCN_mouse_NEMO_snATAC-Seq_combined_cluster_crossSpeciesNeurons_sample_heatmap.pdf', width = 6, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM[c(1,2,3,4,5,6,7,8,10,11,13,14),])), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

###################################################
# 2) plot umap embedding by experiment_short and subclass
pdf('BICCN_mouse_NEMO_snATAC-Seq_combined_umap.pdf', width = 8, height = 5)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "cluster", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
dev.off()