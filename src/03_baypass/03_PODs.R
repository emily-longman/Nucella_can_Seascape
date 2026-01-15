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

# ================================================================================== #

# Create PODs

# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("data/processed/baypass/omega/NC_baypass_summary_beta_params.out", h=T)$Mean

# Omega file
omega <- as.matrix(read.table("data/processed/baypass/omega/NC_baypass_mat_omega.out"))

setwd("data/processed/baypass/PODs")

# Create PODs
foreach(i=1:10, .errorhandling="remove")%do%{
    suffix <- paste("POD.", i, sep="")
    simulate.baypass(omega.mat=omega, nsnp = 500000, beta.pi=pi.beta.coef, sample.size=pooldata@poolsizes, suffix=suffix)

}