suppressMessages(library(ArchR))
options(stringsAsFactors = F)

##################################
### set Arrow File parameters ####
addArchRThreads(threads = 1)
addArchRGenome("hg38")

proj = loadArchRProject("ArchR_BICCN_huMOp_SNARE-Seq2_filtered")

#########################################
# 1) heatmap of cell subclass by experiment_short
library(pheatmap)
cM <- confusionMatrix(paste0(proj$subclass), paste0(proj$patient))

pdf('BICCN_huMOp_SNARE-Seq2_combined_cluster_unnomralized_sample_heatmap.pdf', width = 7, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM)), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

cM <- cM / Matrix::rowSums(cM)

pdf('BICCN_huMOp_SNARE-Seq2_combined_cluster_sample_heatmap.pdf', width = 7, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM)), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

pdf('BICCN_huMOp_SNARE-Seq2_combined_cluster_crossSpeciesNeurons_sample_heatmap.pdf', width = 6, height = 4)
p1 <- pheatmap::pheatmap( mat = as.matrix(t(cM[c(2,3,5,6,7,8,9,10,11,12,15,16,17,18),])), color = paletteContinuous("whiteBlue"), border_color = "black")
dev.off()

###################################################
# 2) plot umap embedding by patient and subclass
pdf('BICCN_huMOp_SNARE-Seq2_combined_umap.pdf', width = 8, height = 5)
p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "patient", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "subclass", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
dev.off()
