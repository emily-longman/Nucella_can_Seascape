# Use poolfstat to convert pooldata to Baypass input files

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

# ================================================================================== #

# Create baypass input files for abiotic mean pH analyses

# Convert to BayPass input file
pooldata2genobaypass(pooldata, writing.dir = "data/processed/baypass/input_files", subsamplesize = -1)
# Three output files = genobaypass (allele counts), poolsize (haploid size per pool), & snpdet (snp info matrix). 
# Subsample size can be used to sample to a smaller number of SNPs. If the subsample size is <0, then all SNPs are included in the BayPass files.

# ================================================================================== #

# Create baypass input files for biotic shell thickness analyses

# Subset data - cut site: 4 (ARA)
pooldata.subset.18pop <- pooldata.subset(pooldata, pool.index=c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset.18pop, writing.dir = "data/processed/baypass/input_files", prefix="subset18pop", subsamplesize = -1)

# ================================================================================== #
