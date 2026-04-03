# Use poolfstat to convert VCF to Baypass input file

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
data_phenotypic_PODs="data/processed/phenotypic_data/PODs"
if (!dir.exists(data_phenotypic_PODs)) {dir.create(data_phenotypic_PODs)}

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Subset data - cutting sites: PSG, KH, FR, PL, PSN, HZD, OCT, STR
pooldata.subset <- pooldata.subset(pooldata, pool.index=c(2,3,4,5,6,7,10,13,15,17,18))

# ================================================================================== #

# Create PODs

# Get estimates (post. mean) of both the a_pi and b_pi parameters of the Pi Beta distribution
pi.beta.coef <- read.table("data/processed/phenotypic_data/omega/NC_pheno_summary_beta_params.out", h=T)$Mean

# Omega file
omega <- as.matrix(read.table("data/processed/phenotypic_data/omega/NC_pheno_mat_omega.out"))

# Set working directory
setwd("data/processed/phenotypic_data/PODs")

# Create PODs - ~8M SNPs (match # SNPs of poolobject)
foreach(i=1:5, .errorhandling="remove")%do%{
    suffix <- paste("POD8M.", i, sep="")
    simulate.baypass(omega.mat=omega, nsnp = pooldata.subset@nsnp, beta.pi=pi.beta.coef, sample.size=pooldata.subset@poolsizes, suffix=suffix)
}
foreach(i=6:10, .errorhandling="remove")%do%{
    suffix <- paste("POD8M.", i, sep="")
    simulate.baypass(omega.mat=omega, nsnp = pooldata.subset@nsnp, beta.pi=pi.beta.coef, sample.size=pooldata.subset@poolsizes, suffix=suffix)
}