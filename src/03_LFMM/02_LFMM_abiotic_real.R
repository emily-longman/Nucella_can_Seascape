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
out_dir_real <- paste("data/processed/LFMM_analyses/real")
if (!dir.exists(out_dir_real)) {dir.create(out_dir_real)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load lfmm data
allele_freqs_t <- fread("data/processed/LFMM_analyses/Ncan_afs.lfmm")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Order bio-oracle data based on pooldata
bio_oracle_sites_2010_ordered <- bio_oracle_sites_2010 %>% mutate(location = factor(location, levels = pooldata@poolnames)) %>% arrange(location)

# Extract just the abiotic variables and scale
abiotic_scaled <- data.frame(scale(bio_oracle_sites_2010_ordered[,4:12]))

# Extract just the variable of interest
abiotic_scaled_var <- abiotic_scaled %>% select(all_of(var))

# ================================================================================== #

# Univariate lfmm

# Run LFMM2 on only one var
mod_lfmm2 <- lfmm2(input = allele_freqs_t, env = abiotic_scaled_var, K = 1, effect.sizes = TRUE)

# Save
save(mod_lfmm2, file = paste0("data/processed/LFMM_analyses/real/mod_lfmm2_", var, ".RData"))

# B - effect size of the environmental variables per SNP
# U - latent factors per pop
# V - loadings of the latent factors on the SNPs

# Compute p-values
pv <- lfmm2.test(object = mod_lfmm2, input = allele_freqs_t, env = abiotic_scaled_var)

# Save
save(pv, file = paste0("data/processed/LFMM_analyses/real/pv_", var,".RData"))

