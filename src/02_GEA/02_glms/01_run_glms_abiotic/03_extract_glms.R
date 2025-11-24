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

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
env_var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_per_env_var/glms_", env_var, sep = "")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis/', pattern = "GLM_100perm_Bio-Oracle_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis/GLM_100perm_Bio-Oracle_chunk_")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(list = ls())
}

