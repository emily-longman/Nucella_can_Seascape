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
pops <- read.table("guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load POD results
XtX.POD.1 <- read.table("data/processed/outlier_analyses/baypass/POD/POD_1/NC_baypass_POD_1_summary_pi_xtx.out", header=T)
pod.1.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/POD_1/NC_baypass_POD_1_summary_beta_params.out", header=T)
XtX.POD.2 <- read.table("data/processed/outlier_analyses/baypass/POD/POD_2/NC_baypass_POD_2_summary_pi_xtx.out", header=T)
pod.2.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/POD_2/NC_baypass_POD_2_summary_beta_params.out", header=T)
XtX.POD.3 <- read.table("data/processed/outlier_analyses/baypass/POD/POD_3/NC_baypass_POD_3_summary_pi_xtx.out", header=T)
pod.3.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/POD_3/NC_baypass_POD_3_summary_beta_params.out", header=T)
XtX.POD.4 <- read.table("data/processed/outlier_analyses/baypass/POD/POD_4/NC_baypass_POD_4_summary_pi_xtx.out", header=T)
pod.4.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/POD_4/NC_baypass_POD_4_summary_beta_params.out", header=T)
XtX.POD.5 <- read.table("data/processed/outlier_analyses/baypass/POD/POD_5/NC_baypass_POD_5_summary_pi_xtx.out", header=T)
pod.5.pi.beta.coef <- read.table("data/processed/outlier_analyses/baypass/POD/POD_5/NC_baypass_POD_5_summary_beta_params.out", header=T)

# ================================================================================== #

# Sanity Check: Compare POD and original data estimates
# Compare the estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the POD analysis

# Extract values of a_beta_pi and b_beta_pi
pi.beta.coef.mean <- pi.beta.coef$Mean
pod.1.pi.beta.coef.mean <- pod.1.pi.beta.coef$Mean
pod.2.pi.beta.coef.mean <- pod.2.pi.beta.coef$Mean
pod.3.pi.beta.coef.mean <- pod.3.pi.beta.coef$Mean
pod.4.pi.beta.coef.mean <- pod.4.pi.beta.coef$Mean
pod.5.pi.beta.coef.mean <- pod.5.pi.beta.coef$Mean

pi.beta.coef.all <- rbind(pi.beta.coef.mean, pod.1.pi.beta.coef.mean, pod.3.pi.beta.coef.mean, pod.4.pi.beta.coef.mean, pod.5.pi.beta.coef.mean)
colnames(pi.beta.coef.all) <- c(" a_pi",  " b_pi")

# Graph pi parameters
pdf("output/figures/outlier_analyses/POD/Baypass_compare_pi_beta_coef_POD.pdf", width = 5, height = 5)
plot(pi.beta.coef.all, xlim=c(0,25), ylim=c(0,25))
abline(a=0,b=1)
dev.off()

# All of the simulations had very similar a_pi and b_pi parameters. 
# The a_pi parameters for the POD are similar to the original estimates, but the b_pi for the POD are very different than the original estimates.

# ================================================================================== #

# Merge baypass results and SNP metadata
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Get the POD XtX

# Compute the 1% threshold (i.e., identify SNPs where the xtx values are above the 99% significance threshold from the POD)
pod.1.thres.99=quantile(XtX.POD.1$M_XtX, probs=0.99)
pod.2.thres.99=quantile(XtX.POD.2$M_XtX, probs=0.99)
pod.3.thres.99=quantile(XtX.POD.3$M_XtX, probs=0.99)
pod.4.thres.99=quantile(XtX.POD.4$M_XtX, probs=0.99)
pod.5.thres.99=quantile(XtX.POD.5$M_XtX, probs=0.99)

# Take the mean of the 5 runs to determine the pod threshold
pod.thres.99 <- mean(c(pod.1.thres.99, pod.2.thres.99, pod.3.thres.99, pod.4.thres.99, pod.5.thres.99))

# Identify outliers -- 318,349 SNPs
baypass_POD_sig_SNPs_threshold_1 <- SNP.XtX.dt[which(SNP.XtX.dt$M_XtX >= pod.thres.99),]

# ================================================================================== #

# Compute the 0.01% threshold (i.e., identify SNPs where the xtx values are above the 99.99% significance threshold from the POD)
pod.1.thres.9999=quantile(XtX.POD.1$M_XtX, probs=0.9999)
pod.2.thres.9999=quantile(XtX.POD.2$M_XtX, probs=0.9999)
pod.3.thres.9999=quantile(XtX.POD.3$M_XtX, probs=0.9999)
pod.4.thres.9999=quantile(XtX.POD.4$M_XtX, probs=0.9999)
pod.5.thres.9999=quantile(XtX.POD.5$M_XtX, probs=0.9999)

# Take the mean of the 5 runs to determine the pod threshold 
pod.thres.9999 <- mean(c(pod.1.thres.9999, pod.2.thres.9999, pod.3.thres.9999, pod.4.thres.9999, pod.5.thres.9999))

# Identify outliers -- 3,095 SNPs
baypass_POD_sig_SNPs_threshold_0.01 <- SNP.XtX.dt[which(SNP.XtX.dt$M_XtX >= pod.thres.9999),]


# ================================================================================== #

# Compute the 0.001% threshold (i.e., identify SNPs where the xtx values are above the 99.999% significance threshold from the POD)
pod.1.thres.99999=quantile(XtX.POD.1$M_XtX, probs=0.99999)
pod.2.thres.99999=quantile(XtX.POD.2$M_XtX, probs=0.99999)
pod.3.thres.99999=quantile(XtX.POD.3$M_XtX, probs=0.99999)
pod.4.thres.99999=quantile(XtX.POD.4$M_XtX, probs=0.99999)
pod.5.thres.99999=quantile(XtX.POD.5$M_XtX, probs=0.99999)

# Take the mean of the 5 runs to determine the pod threshold 
pod.thres.99999 <- mean(c(pod.1.thres.99999, pod.2.thres.99999, pod.3.thres.99999, pod.4.thres.99999, pod.5.thres.99999))

# Identify outliers -- 339 SNPs
baypass_POD_sig_SNPs_threshold_0.001 <- SNP.XtX.dt[which(SNP.XtX.dt$M_XtX >= pod.thres.99999),]


# ================================================================================== #

# Add the threshold to the actual XtX plot
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres_simp.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", pch=19)
abline(h=pod.thres.99, lty=2, col='red') 
dev.off()

# Color by significance
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", col=ifelse(SNP.XtX.dt$M_XtX >= pod.thres.99, "#bebebe93", "black"), pch=19)
abline(h=pod.thres.99, lty=2, col='red') 
dev.off()


# Add the threshold to the actual XtX plot
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres_simp_0.99999.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", pch=19)
abline(h=pod.thres.99999, lty=2, col='red') 
dev.off()

# Color by significance
pdf("output/figures/outlier_analyses/POD/Baypass_xtx_POD_thres_0.99999.pdf", width = 10, height = 5)
plot(SNP.XtX.dt$M_XtX, xlab="Position", col=ifelse(SNP.XtX.dt$M_XtX >= pod.thres.99999, "#bebebe93", "black"), pch=19)
abline(h=pod.thres.99999, lty=2, col='red') 
dev.off()

# ================================================================================== #

# Write file of POD significant SNPs
write.table(baypass_POD_sig_SNPs_threshold_1, "data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_1", sep = "\t", row.names=FALSE)
write.table(baypass_POD_sig_SNPs_threshold_0.01, "data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.01", sep = "\t", row.names=FALSE)
write.table(baypass_POD_sig_SNPs_threshold_0.001, "data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.001", sep = "\t", row.names=FALSE)
