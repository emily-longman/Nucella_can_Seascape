# Graph Baypass seascape - abiotic output.

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen'))
library(data.table)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(geigen)

# ================================================================================== #

# Read in Baypass results
baypass <- read.table("data/processed/GEA/baypass/abiotic/aux_mode/NC_baypass_abiotic_aux_summary_pi_xtx.out", header=T)
betai_snp <- read.table("data/processed/GEA/baypass/abiotic/aux_mode/NC_baypass_abiotic_aux_summary_betai.out", header=T)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #
# ================================================================================== #

# Analyze and graph the baypass results

# Look at the behavior of the p-values associated to the XtX estimator
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_xtx_hist.pdf", width = 5, height = 5)
hist(10**(-1*baypass$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_xtx.pdf", width = 10, height = 5)
plot(baypass$XtXst)
dev.off()

# Graph corrected xtx
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_corrected_xtx.pdf", width = 10, height = 5)
plot(baypass$M_XtX, xlab="SNP", ylab="XtX Corrected")
dev.off()

# Graph outliers
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_xtx_outliers.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) # Adjust margins
plot(baypass$log10.1.pval., ylab="-log10(P-value)", xlab="Position")
abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
abline(h=-log10(0.05/dim(baypass)[1]), lty=1, col="red") # Bonferroni threshold
dev.off()

# ================================================================================== #

# Fancier graphing with ggplot

# Join baypass morphology PC results with snp metadata
baypass_pos <- cbind(snp.meta, baypass)

# Graph with ggplot
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_xtx_outliers_ggplot.pdf", width = 10, height = 5)
ggplot(baypass_pos, aes(x = chr, y = log10.1.pval.)) +
geom_jitter(size = 1, alpha = 0.7, show.legend = FALSE, height=NULL) + 
geom_hline(yintercept = -log10(0.05/dim(baypass)[1]), linetype = "dashed", color = "red", linewidth = 0.8)  +  
theme_classic(base_size = 20) +  
theme(panel.spacing = unit(0.5, "lines"),
axis.text.x = element_blank(), axis.ticks.x = element_blank(),
axis.text.y = element_text(size = 12)) +
labs(x = "Position", y = expression(-log[10](italic(p))))
dev.off()

# ================================================================================== #

# Identify SNP list for baypass results

# Identify significant SNPs (p<0.001)
baypass_pos_sig_SNPs <- data.frame(baypass_pos[which(baypass_pos$log10.1.pval. >= -log10(0.001)),])

# Write file of outlier SNPs
write.table(baypass_pos_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_sig_SNPs.txt", sep="\t", row.names=FALSE)
write.csv(baypass_pos_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_sig_SNPs.csv", row.names=FALSE)
write.csv(baypass_pos_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_sig_SNPs", row.names=FALSE)

# Identify bonferroni significant SNPs
baypass_pos_bonf_sig_SNPs <- data.frame(baypass_pos[which(baypass_pos$log10.1.pval. >= -log10(0.05/dim(baypass_pos)[1])),])

# Write file of bonferroni outlier SNPs
write.table(baypass_pos_bonf_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_bonf_sig_SNPs.txt", sep="\t", row.names=FALSE)
write.csv(baypass_pos_bonf_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_bonf_sig_SNPs.csv", row.names=FALSE)
write.csv(baypass_pos_bonf_sig_SNPs, "data/processed/GEA/baypass/abiotic/aux_mode/baypass_pos_bonf_sig_SNPs", row.names=FALSE)

# ================================================================================== #
# ================================================================================== #

# Plot the estimates of the Bayes Factor (estimates are given in dB unites (i.e., 10xlog10(BF)))
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_estimates_Bayes_Factor.pdf", width = 5, height = 5)
plot(betai_snp$BF.dB., xlab="SNP", ylab="Bayes Factor (in dB)")
dev.off()

# Plot the underlying regression coefficient
pdf("output/figures/GEA/baypass/abiotic/Baypass_aux_regression_coef.pdf", width = 5, height = 5)
plot(betai_snp$M_Beta, xlab="SNP", ylab=expression(beta~"coefficient"))
dev.off()