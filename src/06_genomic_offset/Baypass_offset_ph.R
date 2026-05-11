# Try genomic offset with baypass

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)
library(poolfstat)
library(RColorBrewer)
library(viridis)

# Baypass functions
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")

# ================================================================================== #

# Load data

# Load regression file
regfile <- read.table("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run1_summary_betai_reg.out", header=T)

# Load covariable file
ph.cov.file <- read.table("guide_files/Baypass_ph_mean.txt", header = F)

ph.cov.file.t <- t(ph.cov.file)

ph.future.file <- read.table()

GO <- compute_genetic_offset(regfile = regfile, covfile = ph.cov.file.t)