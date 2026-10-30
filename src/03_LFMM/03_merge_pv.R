# Extract and merge LFMM output

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'magrittr'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(magrittr)

#BiocManager::install("LEA")
library(LEA)


# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)
# Extract specific vars
abiotic_vars <- var_full_names$variable[which(var_full_names$Variable == "Abiotic")]
biotic_MARINe_vars <- var_full_names$variable[which(var_full_names$Variable == "Biotic")][-c(15,16,17,18)]
biotic_Mcali_vars <- var_full_names$variable[which(var_full_names$Variable == "Biotic")][c(15,16,17,18)]

# Path to files
path = paste0('data/processed/LFMM_analyses/real/pv_')

# ================================================================================== #

# Extract p-value for each environmental variable - abiotic

# Snp info
load("data/processed/LFMM_analyses/snp.info.Rdata")
str(snp.info)

# Read all the files and add a column with the chunk
LFMM_real_abiotic <- foreach(w=abiotic_vars, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(paste0(path, w, ".RData")))
    str(o)

    # Extract p-values and zscores per SNP
    tmp <- data.frame(snp.info, o$pvalues, o$zscores)

    # Add column with variable
    tmp$var = var

    # Return output
    tmp
}
# ================================================================================== #

# Extract p-value for each environmental variable - biotic MARINe

# Snp info
load("data/processed/LFMM_analyses/snp.info.16pop.Rdata")
str(snp.info.16pop)

# Read all the files and add a column with the chunk
LFMM_real_biotic_MARINe <- foreach(w=biotic_MARINe_vars, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(paste0(path, w, ".RData")))
    str(o)

    # Extract p-values and zscores per SNP
    tmp <- data.frame(snp.info.16pop, o$pvalues, o$zscores)

    # Add column with variable
    tmp$var = var

    # Return output
    tmp
}

# ================================================================================== #

# Extract p-value for each environmental variable - biotic Mcali

# Snp info
load("data/processed/LFMM_analyses/snp.info.18pop.Rdata")
str(snp.info.18pop)

# Read all the files and add a column with the chunk
LFMM_real_biotic_Mcali <- foreach(w=biotic_Mcali_vars, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(paste0(path, w, ".RData")))
    str(o)

    # Extract p-values and zscores per SNP
    tmp <- data.frame(snp.info.18pop, o$pvalues, o$zscores)

    # Add column with variable
    tmp$var = var

    # Return output
    tmp
}

# ================================================================================== #

# Join all
LFMM_all <- rbind(LFMM_real_abiotic, LFMM_real_biotic_MARINe, LFMM_real_biotic_Mcali)


# Save output
save(LFMM_all, file = "data/processed/LFMM_analyses/LFMM_all.RData")
