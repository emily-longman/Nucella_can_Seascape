# Analyze bbiotic stressors with LFMM

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
out_dir_real <- paste("data/processed/LFMM_analyses/real")
if (!dir.exists(out_dir_real)) {dir.create(out_dir_real)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load lfmm data
allele_freqs_t_16pop <- fread("data/processed/LFMM_analyses/Ncan_afs_16pop.lfmm")

# Load MARINe data
MARINe_data <- read.csv("data/processed/enviro_data/MARINe/MARINe_data_16pops.csv", header=T)

# ================================================================================== #

# Order data based on pooldata
MARINe_data_ordered <- MARINe_data %>% mutate(location = factor(location, levels = pooldata@poolnames)) %>% arrange(location)

# Extract just the abiotic variables and scale
biotic_scaled <- data.frame(scale(MARINe_data_ordered[,13:26]))

# Extract just the variable of interest
biotic_scaled_var <- biotic_scaled %>% select(all_of(var))

# ================================================================================== #

# Identify SNPs with no variance
snp.var <- apply(allele_freqs_t_16pop, 2, stats::var)
# How many SNPs have no variance?
sum(snp.var == 0)
# Index
keep <- snp.var > 0

# Filter SNPs with no variance
allele_freqs_t_16pop_filt <- as.matrix(allele_freqs_t_16pop)[, keep]
# Check dim
dim(allele_freqs_t_16pop_filt)

# ================================================================================== #

# Univariate lfmm

# Do 5 runs
foreach(i=1:5)%do%{

    # Run LFMM2 on only one var
    mod_lfmm2 <- lfmm2(input = allele_freqs_t_16pop_filt, env = biotic_scaled_var, K = 1, effect.sizes = TRUE)

    # Save
    save(mod_lfmm2, file = paste0("data/processed/LFMM_analyses/real/mod_lfmm2_", var, "run_", i, ".RData"))

    # B - effect size of the environmental variables per SNP
    # U - latent factors per pop
    # V - loadings of the latent factors on the SNPs

    # Compute p-values
    pv <- lfmm2.test(object = mod_lfmm2, input = allele_freqs_t_16pop_filt, env = biotic_scaled_var)

    # Save
    save(pv, file = paste0("data/processed/LFMM_analyses/real/pv_", var, "run_", i, ".RData"))

}