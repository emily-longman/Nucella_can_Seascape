# Graph Baypass XtX output with POD.

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'ggplot2'))
library(data.table)
library(ggplot2)

# ================================================================================== #

# Read in Baypass results
XtX <- read.table("data/processed/outlier_analyses/baypass/xtx/NC_baypass_core_summary_pi_xtx.out", header=T)
pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/xtx/NC_baypass_core_summary_beta_params.out", header=T)
XtX.POD <- read.table("data/processed/outlier_analyses/baypass/POD/NC_baypass_POD_summary_pi_xtx.out", header=T)
pod.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/NC_baypass_POD_summary_beta_params.out", header=T)
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #

# Sanity Check: Compare POD and original data estimates

# Compare the estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the POD analysis
pi.beta.coef.mean <- pi.beta.coef$Mean
pod.pi.beta.coef.mean <- pod.pi.beta.coef$Mean

# Graph pi parameters
pdf("output/figures/outlier_analyses/POD/Baypass_compare_pi_beta_coef_POD.pdf", width = 5, height = 5)
plot(pod.pi.beta.coef.mean, pi.beta.coef.mean)
abline(a=0,b=1)
dev.off()

# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Get the POD XtX

# Compute the 5% threshold (i.e., identify SNPs where the xtx values are above the 99% significance threshold from the POD)
pod.thres.99=quantile(XtX.POD$M_XtX, probs=0.99)

# Identify outliers -- 320,742 SNPs
baypass_POD_sig_SNPs <- SNP.XtX.dt[which(SNP.XtX.dt$M_XtX >= pod.thres.99),]

# Add the threshold to the actual XtX plot
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres_simp.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", pch=19)
abline(h=pod.thres.99, lty=2, col='red') 
dev.off()

# Color by significance
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", col=ifelse(SNP.XtX.dt$M_XtX >= pod.thres.99, "black", "#bebebe93"), pch=19)
abline(h=pod.thres.99, lty=2, col='red') 
dev.off()

# More detailed graphing -- this doesn't look correct
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_outliers.pdf", width = 10, height = 5)
par(mar=c(5,5,4,1)+.1) #Adjust margins
plot(SNP.XtX.dt$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position", 
col=ifelse(SNP.XtX.dt$M_XtX >= pod.thres.99, "black", "#bebebe93"), pch=19)
abline(h=pod.thres.99, lty=2, col='red') 
abline(h=-log10(0.05/dim(SNP.XtX.dt)[1]), lty=1, col="blue") # Bonferroni threshold
dev.off()

# ================================================================================== #

# Write file of POD significant SNPs
write.table(baypass_POD_sig_SNPs, "data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs", sep = "\t", row.names=FALSE)
