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
data_processed_outlier="data/processed/baypass/PODs"
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

# Create PODs

# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("data/processed/baypass/omega/NC_baypass_summary_beta_params.out", h=T)$Mean

# Omega file
omega <- as.matrix(read.table("data/processed/baypass/omega/NC_baypass_mat_omega.out"))

# Create PODs
POD.sim <- simulate.baypass(omega.mat=omega, nsnp = pooldata@nsnp, beta.pi=pi.beta.coef, sample.size=pooldata@poolsizes, suffix="NC_POD_")