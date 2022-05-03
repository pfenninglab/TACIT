# Motor cortex, solitary
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_solitary_results_perm1m.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.005222
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_solitary_results_perm1m_adjp.csv", row.names = FALSE)

# Motor cortex, solitary, genome-wide
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_full_solitary_results_perm1m_filt.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.5598414
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_full_solitary_results_perm1m_filt_adjp.csv", row.names = FALSE, quote = FALSE)

# Motor cortex, group living
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_groupliving_results_perm100k.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.00826
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_groupliving_results_perm100k_adjp.csv", row.names = FALSE)

# Motor cortex, group living, genome-wide
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/cortex_full_groupLiving_results_perm1m_filt.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.20052
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/cortex_full_groupLiving_results_perm1m_filt_adjp.csv", row.names = FALSE, quote = FALSE)

# Liver, solitary
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/liver_solitary_results_perm10k.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.4267
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/liver_solitary_results_perm10k_adjp.csv", row.names = FALSE, quote = FALSE)

# Liver, group living
results = read.csv(file="/projects/MPRA/Irene/PGLSOutputs/liver_groupLiving_results_perm10k.csv")
rownames(results) = results$Enhancer
bh = p.adjust(results$Exp_Pvalue, method = "BH")
min(bh) # 0.18837
results$BH = bh
write.csv(results, "/projects/MPRA/Irene/PGLSOutputs/liver_groupLiving_results_perm10k_adjp.csv", row.names = FALSE, quote = FALSE)