# Graph Baypass C2 output.

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
NC.S.vs.N.contrast.C2 <- read.table("data/processed/outlier_analyses/baypass/C2/NC_baypass_C2_summary_contrast.out", header=T)

NC.S.vs.N.contrast.C2.no_PGP <- read.table("data/processed/outlier_analyses/baypass/C2/NC_baypass_C2_no_PGP_summary_contrast.out", header=T)

# Note: on baypass manual, they also looked at "outprefix_summary_betai_reg.out", which contains estaimtes of the Bayes Factor associated with each SNP for each population covariable


# ================================================================================== #


# Check the behavior of the p-values associated to the C2
pdf("output/figures/outlier_analyses/Baypass_C2_hist.pdf", width = 5, height = 5)
hist(10**(-1*NC.S.vs.N.contrast.C2$log10.1.pval.),freq=F,breaks=50)
abline(h=1, col="red")
dev.off()

pdf("output/figures/outlier_analyses/Baypass_C2_outliers.pdf", width = 5, height = 5)
plot(NC.S.vs.N.contrast.C2$log10.1.pval., ylab="C2 p-value (-log10 scale)")
abline(h=3,lty=2) #0.001 p--value theshold
dev.off()

# Example from baypass manual pg 35 - comparing Bayes Factor with C2 results
#plot(lsa.ecotype.bf,lsa.ecotype.C2$log10.1.pval.,
#xlab="BF",ylab="C2 p-value (-log10 scale)")
#abline(h=3,lty=2) #0.001 p--value theshold
#abline(v=20,lty=2) #BF threshold for decisive evidence (according to Jeffreys’ rule)


# No PGP in contrast

# Check the behavior of the p-values associated to the C2
pdf("output/figures/outlier_analyses/Baypass_C2_hist_no_PGP.pdf", width = 5, height = 5)
hist(10**(-1*NC.S.vs.N.contrast.C2.no_PGP$log10.1.pval.),freq=F,breaks=50)
abline(h=1, col="red")
dev.off()

pdf("output/figures/outlier_analyses/Baypass_C2_outliers_no_PGP.pdf", width = 5, height = 5)
plot(NC.S.vs.N.contrast.C2.no_PGP$log10.1.pval., ylab="C2 p-value (-log10 scale)")
abline(h=3,lty=2, col="red") #0.001 p--value theshold
dev.off()
