# Merge glms output

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
#install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/real/', pattern = "GLM_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/real/"), x))))

# Check number of files
length(file_names_v)
file_names_v=file_names_v[1:3]
# Read all the files and add a column with the chunk
glm.model.collated.real =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.model.collated.real)

# ================================================================================== #

# Save merged data
save(glm.model.collated.real, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.real.Rdata"))
