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
out_dir_perm <- paste("data/processed/LFMM_analyses/perm")
if (!dir.exists(out_dir_perm)) {dir.create(out_dir_perm)}
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

# Extract and scale variable of interest
scaled_var <- MARINe_data %>% select(all_of(var)) %>% scale()

# Permutations to generate null expectation for LFMM
foreach(l=61:80)%do%{

    # Set seed
    set.seed(l)
        
    # Shuffle enviro data for 'j' enviro variable
    scaled_var_shuffle <- sample(scaled_var)

    # Univariate lfmm

    # Run LFMM2 on only one var
    mod_lfmm2 <- lfmm2(input = allele_freqs_t_16pop_filt, env = scaled_var_shuffle, K = 1, effect.sizes = TRUE)

    # Save
    save(mod_lfmm2, file = paste0("data/processed/LFMM_analyses/perm/mod_lfmm2_", var, "_perm", l, ".RData"))

    # Compute p-values
    pv <- lfmm2.test(object = mod_lfmm2, input = allele_freqs_t_16pop_filt, env = scaled_var_shuffle)

    # Save
    save(pv, file = paste0("data/processed/LFMM_analyses/perm/pv_", var,"_perm", l, ".RData"))

    # Housekeeping
    rm(mod_lfmm2)
    rm(pv)

}
