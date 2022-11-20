
pvals = read.csv("brainsize_pdists.csv")

svg(file = "cortex_p.svg", width=5, height=5)
d = density(pvals$Cortex_Perm_Pvalue)
plot(d, xlab = "Permulations p-values", ylab="Density", main="",ylim=c(0,2))
polygon(d, col="grey", border="black", lwd=2)
dev.off()

svg(file = "pv_p.svg", width=5, height=5)
d = density(na.omit(pvals$PV_Perm_Pvalue))
plot(d, xlab = "Permulations p-values", ylab="Density", main="", ylim=c(0,2))
polygon(d, col="grey", border="black", lwd=2)
dev.off()