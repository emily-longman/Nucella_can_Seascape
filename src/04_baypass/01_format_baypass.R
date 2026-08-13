# Use poolfstat to convert VCF to Baypass input file and create PODs

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
install.packages(c('poolfstat'))
library(poolfstat)

# Baypass functions
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")

# ================================================================================== #

# Generate Folders and files

# Make output directories
data_processed_outlier="data/processed/baypass"
if (!dir.exists(data_processed_outlier)) {dir.create(data_processed_outlier)}
data_processed_outlier="data/processed/baypass/input_files"
if (!dir.exists(data_processed_outlier)) {dir.create(data_processed_outlier)}

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Convert to BayPass input file
pooldata2genobaypass(pooldata, writing.dir = "data/processed/baypass/input_files", subsamplesize = -1)
# Three output files = genobaypass (allele counts), poolsize (haploid size per pool), & snpdet (snp info matrix). 
# Subsample size can be used to sample to a smaller number of SNPs. If the subsample size is <0, then all SNPs are included in the BayPass files.

# ================================================================================== #

# Subset data - cutting sites: 6 (CBL), 15 (VD), 16 (OCT)
pooldata.subset <- pooldata.subset(pooldata, pool.index=c(1,2,3,4,5,7,8,9,10,11,12,13,14,17,18,19))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset, writing.dir = "data/processed/baypass/input_files", prefix="subset", subsamplesize = -1)

# ================================================================================== #

# Subset data - cut site: 4 (ARA)
pooldata.subset.18pop <- pooldata.subset(pooldata, pool.index=c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset.18pop, writing.dir = "data/processed/baypass/input_files", prefix="subset18pop", subsamplesize = -1)

# ================================================================================== #


# Subset data - cut site: 4 (ARA) and 9 (HZD)
pooldata.subset.17pop <- pooldata.subset(pooldata, pool.index=c(1,2,3,5,6,7,8,10,11,12,13,14,15,16,17,18,19))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset.17pop, writing.dir = "data/processed/baypass/input_files", prefix="subset17pop", subsamplesize = -1)

# ================================================================================== #


# 9 pop Subset data - cut site: 6 (CBL), 15 (VD), 16 (OCT) AND 1 (STR), 9 (HZD), 10 (SBR), 14 (PL), 17 (PB), 18 (PGP), 19 (PSN)
pooldata.subset.9pop <- pooldata.subset(pooldata, pool.index=c(2,3,4,5,7,8,11,12,13))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset.9pop, writing.dir = "data/processed/baypass/input_files", prefix="subset9pop", subsamplesize = -1)
