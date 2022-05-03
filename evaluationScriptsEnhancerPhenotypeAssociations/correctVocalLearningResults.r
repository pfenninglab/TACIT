results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.03914269
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_adjp.csv", row.names = FALSE)