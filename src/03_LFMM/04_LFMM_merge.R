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



# Extract p-value for each environmental variable - real

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/LFMM_analyses/real/', pattern = "pv_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/LFMM_analyses/real/', x))))

# Extract just abiotic
file_names_v <- file_names_v[c(5,18,21,22,23,24,25,26,27)]

# Read all the files and add a column with the chunk
LFMM_real <- foreach(w=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Extract p-values and zscores per SNP
    tmp <- data.frame(o$pvalues, o$zscores)

    # Add column with variable
    tmp %<>% mutate(var = w) %>% mutate(var = str_remove(var, pattern = "data/processed/LFMM_analyses/real/pv_*"))

    # Remove end of var name
    tmp %<>% mutate(var = str_remove(var, pattern = ".RData"))

    # Create columns specifying perm versus real
    tmp$data = "real"
    tmp$perm = 0

    # Return output
    tmp
}

# Save output
#save(LFMM_real, file = "data/processed/LFMM_analyses/LFMM_real_all.RData")
save(LFMM_real, file = "data/processed/LFMM_analyses/LFMM_real_abiotic.RData")

