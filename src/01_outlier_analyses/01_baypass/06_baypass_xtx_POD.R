# Use pseudo-observed data (POD) to calibrate the XtX

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify path from root
path_from_root <- find_root_file("data", "processed", "outlier_analyses", "baypass", criterion = has_file("README.md"))
# Set working directory as path from root
setwd(path_from_root)

# ================================================================================== #

# Load packages
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages('mvtnorm')
library(mvtnorm)

# ================================================================================== #

# Read in Baypass results

# Read in baypass input gfile 
NC.genobaypass <- geno2YN("genobaypass")

# Read in Baypass results
NC.omega <- as.matrix(read.table("omega/NC_baypass_mat_omega.out"))
beta_params <- read.table("xtx/NC_baypass_core_summary_beta_params.out", header=T)

# ================================================================================== #

# Extract posterior means of a_pi and b_pi
beta_params_mean <- beta_params$Mean

# Create the POD
sim.NC <- simulate.baypass(omega.mat=NC.omega, nsnp=8277206, sample.size=NC.genobaypass$NN, beta.pi=beta_params_mean, pi.maf=0, suffix="NC.baypass.sim")

# NOTE: the above code produces a simulated gfile (G.NC.baypass.sim) in the current folder.