results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m.csv")
rownames(results) = results$Enhancer
results$Orth_count=0
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.03914269
results$bh = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_adjp.csv", row.names = FALSE)

# Use only species with at least 16 orthologs in vocal learners
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_atLeast16True.csv")
rownames(results) = results$Enhancer
results$Orth_count=0
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.0356712
results$bh = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_adjp_atLeast16True.csv", row.names = FALSE)

# Use only species with at least 3 orthologs in vocal learners and orthologs in vocal learners from at least 2 orders
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.03761615
results$bh = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_adjp_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner.csv", row.names = FALSE)

# Use only species with at least 3 orthologs in vocal learners and orthologs in vocal learners from at least 2 orders with repeated peaks removed
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner_noRep.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.03632462
results$bh = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_VocalLearner_results_perm1m_adjp_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner_noRep.csv", row.names = FALSE)

# Use only species with at least 3 orthologs in vocal learners and orthologs in vocal learners from at least 2 orders with repeated peaks removed and corrected p-values
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_vl_results_combined_mousefix_12-31_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.023611
results$bh = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_vl_results_combined_mousefix_12-31_atLeast3True_orthologInChiropteraVocalLearner_orthologInNonChiropteraVocalLearner_adjp.csv", row.names = FALSE)