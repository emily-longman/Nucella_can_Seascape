# Create windows

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

# Load and merge data

# Create list of file names
path <- paste("/data/processed/baypass/window_summary/window_chunk_analysis_Mtross_mean")
file_names = as.list(dir(path = path, pattern = "window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("/gpfs3/scratch/elongman/glms_per_env_var/glms_", env_var, "/", sep = ""), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
glm.model.collated =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.model.collated)

# ================================================================================== #

# Save merged data
save(glm.model.collated, file = paste("/gpfs3/scratch/elongman/glms_per_env_var/glm.collated_", env_var,".Rdata", sep = ""))
