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
data_phenotypic_data="data/processed/phenotypic_data"
if (!dir.exists(data_phenotypic_data)) {dir.create(data_phenotypic_data)}
data_phenotypic_data_baypass="data/processed/phenotypic_data/baypass_input_files"
if (!dir.exists(data_phenotypic_data_baypass)) {dir.create(data_phenotypic_data_baypass)}

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Subset data - cutting sites: PSG, KH, FR, PL, PSN, HZD, OCT, STR
pooldata.subset <- pooldata.subset(pooldata, pool.index=c(2,3,4,5,6,7,10,13,15,17,18))

# Convert to BayPass input file
pooldata2genobaypass(pooldata.subset, writing.dir = "data/processed/phenotypic_data/baypass_input_files", subsamplesize = -1)
# Three output files = genobaypass (allele counts), poolsize (haploid size per pool), & snpdet (snp info matrix). 
# Subsample size can be used to sample to a smaller number of SNPs. If the subsample size is <0, then all SNPs are included in the BayPass files.

# ================================================================================== #
