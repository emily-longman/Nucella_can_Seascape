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

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Extract p-value for each environmental variable - real

# Snp info
load("data/processed/LFMM_analyses/snp.info.16pop.Rdata")
str(snp.info.16pop)

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/LFMM_analyses/real/', pattern = paste0("pv_", var, "run_*")))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/LFMM_analyses/real/', x))))

# Read all the files and add a column with the chunk
LFMM_real <- foreach(w=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))
    str(o)

    # Extract p-values and zscores per SNP
    tmp <- data.frame(snp.info.16pop, o$pvalues, o$zscores)

    # Add column with variable
    tmp$var = var

    # Return output
    tmp
}

# ================================================================================== #

# Save output
save(LFMM_real, file = paste0("data/processed/LFMM_analyses/LFMM_", var, ".RData"))
