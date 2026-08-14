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
install.packages(c('poolfstat', 'foreach'))
library(poolfstat)
library(foreach)

# Baypass functions
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")

# ================================================================================== #

# Generate Folders and files

# Make output directories
data_processed_outlier="data/processed/baypass/PODs"
if (!dir.exists(data_processed_outlier)) {dir.create(data_processed_outlier)}

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Subset pooldata
pooldata.subset <- pooldata.subset(pooldata, pool.index=c(1,2,3,4,5,7,8,9,10,11,12,13,14,17,18,19))

# ================================================================================== #

# Create PODs

# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("data/processed/baypass/omega_subset/NC_subset_baypass_summary_beta_params.out", h=T)$Mean

# Omega file
omega <- as.matrix(read.table("data/processed/baypass/omega_subset/NC_subset_baypass_mat_omega.out"))

# Set working directory
setwd("data/processed/baypass/PODs")

# Create PODs - ~8M SNPs (match # SNPs of poolobject)
foreach(i=1:5, .errorhandling="remove")%do%{
    suffix <- paste("POD8M_subset.", i, sep="")
    simulate.baypass(omega.mat=omega, nsnp = pooldata.subset@nsnp, beta.pi=pi.beta.coef, sample.size=pooldata.subset@poolsizes, suffix=suffix)
}
foreach(i=6:10, .errorhandling="remove")%do%{
    suffix <- paste("POD8M_subset.", i, sep="")
    simulate.baypass(omega.mat=omega, nsnp = pooldata.subset@nsnp, beta.pi=pi.beta.coef, sample.size=pooldata.subset@poolsizes, suffix=suffix)
}