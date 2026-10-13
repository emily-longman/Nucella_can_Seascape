# Analyze abiotic stressors with LFMM

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

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_dir_biotic <- paste("data/processed/LFMM_analyses/biotic")
if (!dir.exists(out_dir_biotic)) {dir.create(out_dir_biotic)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load lfmm data
allele_freqs_t_18pop <- fread("data/processed/LFMM_analyses/Ncan_afs_18pop.lfmm")

# Load M californianus shell thickness data
Mcalifornianus_data <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean_18pop.csv", header=T)

# ================================================================================== #

# Order data based on pooldata
Mcalifornianus_data_ordered <- Mcalifornianus_data %>% mutate(Site.Code = factor(Site.Code, levels = pooldata@poolnames)) %>% arrange(Site.Code)

# Extract just the abiotic variables and scale
biotic_scaled <- data.frame(scale(Mcalifornianus_data_ordered[,6:9]))

# Extract just the variable of interest
biotic_scaled_var <- biotic_scaled %>% select(all_of(var))

# ================================================================================== #

# Univariate lfmm

# Run LFMM2 on only ph mean
mod_lfmm2 <- lfmm2(input = allele_freqs_t_18pop, env = biotic_scaled_var, K = 1, effect.sizes = TRUE)

# Save
save(mod_lfmm2, file = paste0("data/processed/LFMM_analyses/biotic/mod_lfmm2_", var, ".RData"))

# B - effect size of the environmental variables per SNP
# U - latent factors per pop
# V - loadings of the latent factors on the SNPs

# Compute p-values
pv <- lfmm2.test(object = mod_lfmm2, input = allele_freqs_t_18pop, env = biotic_scaled_var)

# Save
save(pv, file = paste0("data/processed/LFMM_analyses/abiotic/pv_", var,".RData"))

