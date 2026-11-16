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

# Extract p-value for each environmental variable - perm

# Snp info
load("data/processed/LFMM_analyses/snp.info.Rdata")
str(snp.info)

# Create list of file names for perm data
file_names = as.list(dir(path = 'data/processed/LFMM_analyses/perm/', pattern = paste0("pv_", var, "_*")))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/LFMM_analyses/perm/', x))))

# Read all the files and add a column with the chunk
LFMM_perm <- foreach(w=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Extract p-values and zscores per SNP
    tmp <- data.frame(snp.info, o$pvalues, o$zscores)

    # Add column with variable
    tmp$var = var

    # Add column with permutation and perm number
    tmp$data = "permutation"
    tmp %<>% mutate(perm = w) %>% mutate(perm = str_remove(perm, pattern = paste0("data/processed/LFMM_analyses/perm/pv_", var,"_perm"))) 
    # Remove end of var name
    tmp %<>% mutate(perm = str_remove(perm, pattern = ".RData"))

    # Return output
    tmp
}

# Save output
save(LFMM_perm, file = paste0("data/processed/LFMM_analyses/LFMM_perm_", var, ".RData"))
