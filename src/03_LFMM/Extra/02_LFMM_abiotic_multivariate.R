# Analyze abiotic stressors with LFMM

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
#install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)

#BiocManager::install("LEA")
library(LEA)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load lfmm data
allele_freqs_t <- fread("data/processed/LFMM_analyses/Ncan_afs.lfmm")

# Load SNP info
load("data/processed/LFMM_analyses/snp.info.Rdata")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Order bio-oracle data based on pooldata
bio_oracle_sites_2010_ordered <- bio_oracle_sites_2010 %>% mutate(location = factor(location, levels = pooldata@poolnames)) %>% arrange(location)

# Extract just the abiotic variables and scale
abiotic_scaled <- data.frame(scale(bio_oracle_sites_2010_ordered[,4:12]))

# ================================================================================== #

# Univariate lfmm - 

# Run LFMM2 on only ph mean
mod_lfmm2_ph_mean <- lfmm2(input = allele_freqs_t, env = abiotic_scaled$ph_mean, K = 1, effect.sizes = TRUE)

# B - effect size of the environmental variables per SNP
# U - latent factors per pop
# V - loadings of the latent factors on the SNPs

# Compute p-values
pv_ph_mean <- lfmm2.test(object = mod_lfmm2_ph_mean, input = allele_freqs_t, env = abiotic_scaled$ph_mean)

# Show p-value distribution
pdf("output/figures/LFMM_analyses/lfmm_pvalue_dist.pdf", width = 12, height = 6)
ggplot(pv_ph_mean, aes(x = pv.pvalues)) + geom_histogram() + theme_linedraw(base_size = 30)
dev.off()

# Show p-value distribution
pdf("output/figures/LFMM_analyses/lfmm_pvalue.pdf", width = 12, height = 6)
ggplot(pv_ph_mean, aes(x = seq_len(nrow(pvalues)), y = -log10(pv.pvalues))) + geom_point(alpha = 0.7) + theme_linedraw(base_size = 30)
dev.off()

# ================================================================================== #

# Merge with SNP info and do window analysis



# ================================================================================== #

# Shuffle




#########

# Multivariate lfmm - including all abiotic variables - however, many variables are correlated so not the best approach

# Run LFMM2 on only ph mean and the test
mod.lfmm2 <- lfmm2(input = allele_freqs_t, env = abiotic_scaled, K = 1, effect.sizes = TRUE)

# Compute p-values
pv <- lfmm2.test(object = mod.lfmm2, input = test, env = abiotic_scaled, full = TRUE)
# Extract p-values
pvalues = pv$pvalues

# Show p-value distribution
pdf("data/processed/GEA/LEA_analyses/lfmm_pvalue_dist.pdf", width = 12, height = 6)
ggplot(pvalues, aes(x = pv.pvalues)) + geom_histogram() + theme_linedraw(base_size = 30)
dev.off()

# Show p-value distribution
pdf("data/processed/GEA/LEA_analyses/lfmm_pvalue.pdf", width = 12, height = 6)
ggplot(pvalues, aes(x = seq_len(nrow(pvalues)), y = -log10(pv.pvalues))) + geom_point(alpha = 0.7) + theme_linedraw(base_size = 30)
dev.off()

# Note: they recommend using several runs and to combine p-value

# Adjust p-values for multiple testing using Benjamini-Hochberg procedure
for (alpha in c(0.01, 0.05, 0.10)){
    print(paste("Expected FDR:", alpha))
    L = length(pvalues)
    # Return a list of candidates with expected FDR alpha
    # Benjamini-Hochberg's algorithm:
    w = which(sort(pvalues) < alpha * (1:L) / L)
    candidates = order(pvalues)[w]
    # estimated FDR and True Positive Rate
    Lc = length(candidates)
    estimated.FDR = sum(candidates <= 350)/Lc
    print(paste("Observed FDR:",
    round(estimated.FDR, digits = 2)))
    estimated.TPR = sum(candidates > 350)/50
    print(paste("Estimated TPR:",
    round(estimated.TPR, digits = 2)))
}


####### 

# Run example in tutorial to see if can extract B for each var
# load simulated data
data("offset_example")
# 200 diploid individuals genotyped at 510 SNP
Y <- offset_example$geno
# 4 environmental variables
X <- offset_example$env
# Run LFMM2 with 2 latent factors
mod.lfmm2.test <- lfmm2(input = Y, env = X, K = 2, effect.sizes = TRUE)
# Extract pvalues
pv <- lfmm2.test(object = mod.lfmm2.test, input = Y, env = X, full = TRUE)