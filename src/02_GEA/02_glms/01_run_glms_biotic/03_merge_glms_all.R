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

out_dir <- paste("data/processed/GEA/glms/glms_output_marine")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/', pattern = "glm.model.collated.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/', x))))

# Read all the files and add a column with the chunk
glm.model.collated =  
foreach(i=file_names_v, .combine="rbind")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
} 

# ================================================================================== #

# Save merged data
save(glm.model.collated, file = "data/processed/GEA/glms/glms_output_marine/glm.model.collated.marine.Rdata")