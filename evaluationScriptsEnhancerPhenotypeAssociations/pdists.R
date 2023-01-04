
pvals = read.csv("brainsize_pdists.csv")

svg(file = "cortex_staggered_p.svg", width=4.5, height=5)
hist(na.omit(pvals$Cortex_BrainSize_Staggered), xlab = "Brain Size Staggered Permulations p-value", ylab="# of Cortex OCRs", main="")
dev.off()

svg(file = "cortex_staggered_p_zoom.svg", width=4.5, height=5)
hist(na.omit(pvals$Cortex_BrainSize_Staggered), xlab = "Brain Size Staggered Permulations p-value", ylab="# of Cortex OCRs", main="", breaks=2000, xlim=c(0,0.01))
dev.off()

svg(file = "cortex_10k_p.svg", width=4.5, height=5)
hist(na.omit(pvals$Cortex_BrainSize_10K), xlab = "Brain Size 10K Permulations p-value", ylab="# of Cortex OCRs", main="")
dev.off()

svg(file = "pv_staggered_p.svg", width=4.5, height=5)
hist(na.omit(pvals$PV_BrainSize_Staggered), xlab = "Brain Size Staggered Permulations p-value", ylab="# of PV Neuron OCRs", main="")
dev.off()

svg(file = "pv_10k_sol_p.svg", width=4.5, height=5)
hist(na.omit(pvals$PV_Solitary_10K), xlab = "Solitary 10K Permulations p-value", ylab="# of PV Neuron OCRs", main="", xlim=c(0,1))
dev.off()
