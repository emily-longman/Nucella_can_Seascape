# Compare abiotic and biotic stressors with LFMM

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

BiocManager::install("LEA")
library(LEA)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/LEA_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load AF per SNP - generated in src/04_baypass/06_AFs_all.R
load("data/processed/outlier_analyses/afs.all.RData")
str(afs.all)

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Format data for LFMM

# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Extract and manipulate count and coverage

# Extract read count data for SNPs
ref_count <- pooldata@refallele.readcount
ref_count %>% as.data.frame -> count
names(count) = c(pooldata@poolnames)
row.names(count) <- snp.info$SNP_id

# Extract coverage data for SNPs
coverage <- pooldata@readcoverage
coverage %>% as.data.frame -> cov
names(cov) = c(pooldata@poolnames)
row.names(cov) <- snp.info$SNP_id

# Calculate allele frequency for SNPs and transpose
allele_freqs <- count/cov

# Transpose
allele_freqs_t <- t(allele_freqs)

# Write
write.lfmm(allele_freqs_t, "data/processed/GEA/LEA_analyses/Ncan_afs.lfmm")

# Subset to only 500 SNPs to test LFMM
test <- allele_freqs_t[,1:500]
write.lfmm(test, "data/processed/GEA/LEA_analyses/Ncan_test.lfmm")

# ================================================================================== #

# Order bio-oracle data based on pooldata
bio_oracle_sites_2010_ordered <- bio_oracle_sites_2010 %>% mutate(location = factor(location, levels = pooldata@poolnames)) %>% arrange(location)

# Extract just ph and scale
ph_mean_scaled <- scale(bio_oracle_sites_2010_ordered$ph_mean)

# Extract just the abiotic variables and scale
abiotic_scaled <- data.frame(scale(bio_oracle_sites_2010_ordered[,4:12]))

# ================================================================================== #

# PCA with LEA
pc = pca("data/processed/GEA/LEA_analyses/Ncan_afs.lfmm", scale = T)

tw = tracy.widom(pc)

# Plot PC1 and PC2 with ggplot
pdf("data/processed/GEA/LEA_analyses/screeplot.pdf", width = 6, height = 6)
ggplot(tw, aes(x=N, y=percentage)) + 
geom_point(size=5) + theme_linedraw(base_size = 30)
dev.off()

# To be consistent with the GLM approach, use 1 latent factor
# Note: "In principle, the latent factors include levels of population structure due to shared demographic history 
# of background genetic variation but they are generally different of population structure estimates (such as principal components)"

# ================================================================================== #

# Univariate lfmm - mean pH


# Run LFMM2 on only ph mean
mod_lfmm2_ph_mean <- lfmm2(input = allele_freqs_t, env = abiotic_scaled$ph_mean, K = 1)

# Compute p-values
pv_ph_mean <- lfmm2.test(object = mod_lfmm2_ph_mean, input = allele_freqs_t, env = abiotic_scaled$ph_mean)
# Extract p-values
pvalues_ph_mean = pv_ph_mean$pvalues

# Show p-value distribution
pdf("data/processed/GEA/LEA_analyses/lfmm_pvalue_dist.pdf", width = 12, height = 6)
ggplot(pvalues, aes(x = pv.pvalues)) + geom_histogram() + theme_linedraw(base_size = 30)
dev.off()

# Show p-value distribution
pdf("data/processed/GEA/LEA_analyses/lfmm_pvalue.pdf", width = 12, height = 6)
ggplot(pvalues, aes(x = seq_len(nrow(pvalues)), y = -log10(pv.pvalues))) + geom_point(alpha = 0.7) + theme_linedraw(base_size = 30)
dev.off()



# Multivariate lfmm - including all abiotic variables

# Run LFMM2 on only ph mean and the test
mod.lfmm2 <- lfmm2(input = test, env = abiotic_scaled, K = 1)

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
mod.lfmm2 <- lfmm2(input = Y, env = X, K = 2)
# Extract pvalues
pv <- lfmm2.test(object = mod.lfmm2, input = Y, env = X, full = TRUE)