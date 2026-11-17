# Format data for LFMM analyses

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

# ================================================================================== #

# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Write snp info
save(snp.info, file="data/processed/LFMM_analyses/snp.info.Rdata")

# ================================================================================== #

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

# Transpose - 
allele_freqs_t <- t(allele_freqs)

# Write
write.lfmm(allele_freqs_t, "data/processed/LFMM_analyses/Ncan_afs.lfmm")

# ================================================================================== #

# Subset for M californianus analyses - i.e., remove ARA population
allele_freqs_t_18pop <- allele_freqs_t[-which(pooldata@poolnames == "ARA"),]

# Write
write.lfmm(allele_freqs_t_18pop, "data/processed/LFMM_analyses/Ncan_afs_18pop.lfmm")

# Identify SNPs with no variance
snp.var <- apply(allele_freqs_t_18pop, 2, stats::var)
# How many SNPs have no variance?
sum(snp.var == 0)
# Index
keep <- snp.var > 0

# Make snp_id column
snp.info.18pop <- snp.info[keep,]

# Write snp info - 8,274,302 SNPs
save(snp.info.18pop, file="data/processed/LFMM_analyses/snp.info.18pop.Rdata")

# ================================================================================== #

# Subset for MARINe analyses - i.e., remove CBL, VD, OCT populations
allele_freqs_t_16pop <- allele_freqs_t[-which(pooldata@poolnames == "CBL" | pooldata@poolnames == "VD" | pooldata@poolnames == "OCT"),]

# Write
write.lfmm(allele_freqs_t_16pop, "data/processed/LFMM_analyses/Ncan_afs_16pop.lfmm")

# Identify SNPs with no variance
snp.var <- apply(allele_freqs_t_16pop, 2, stats::var)
# How many SNPs have no variance?  SNPs
sum(snp.var == 0)
# Index
keep <- snp.var > 0

# Make snp_id column
snp.info.16pop <- snp.info[keep,]

# Write snp info - 8,260,461 SNPs
save(snp.info.16pop, file="data/processed/LFMM_analyses/snp.info.16pop.Rdata")

# ================================================================================== #

# PCA with LEA
pc = pca("data/processed/GEA/LEA_analyses/Ncan_afs.lfmm", scale = T)

tw = tracy.widom(pc)

# Plot Screeplot
pdf("output/figures/LFMM_analyses/screeplot.pdf", width = 6, height = 6)
ggplot(tw, aes(x=N, y=percentage)) + 
geom_point(size=5) + theme_linedraw(base_size = 30)
dev.off()

# To be consistent with the GLM approach, use 1 latent factor
# Note: "In principle, the latent factors include levels of population structure due to shared demographic history 
# of background genetic variation but they are generally different of population structure estimates (such as principal components)"
