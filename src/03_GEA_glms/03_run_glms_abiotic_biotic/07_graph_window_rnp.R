# Graph window rnp

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'doMC'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic/")
file_names = as.list(dir(path = path, pattern = "glm_window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
win.out =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(win.out)

# ================================================================================== #

# Read in SNP data from Baypass
snpdet <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snpdet) <- c("chr", "pos", "allele1", "allele2")
# Make unique list of chr names
snpdet.chr <- unique(snpdet$chr)

# ================================================================================== #

# Make sure windows are ordered in same chr list as snpdet
win.out.order <- win.out[order(factor(win.out$chr, levels = snpdet.chr)),]

# Note: number of chr between win.out.order and snpdet don't match becuase several chr failed the filter when generating the windows

# ================================================================================== #

# Save merged data
save(win.out.order, file = "data/processed/baypass/window_summary/window_analysis_ph_mean.RData")
load("data/processed/baypass/window_summary/window_analysis_ph_mean.RData")
