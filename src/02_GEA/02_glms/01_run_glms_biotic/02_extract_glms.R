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

# Extract glm chunk for each environmental variable - real

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/real/', pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/real/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_marine/real/GLM_MARINe_chunk_*")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_real", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(glm.model.output)
}


# ================================================================================== #

# Extract glm chunk for each environmental variable - perm 1:25

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/perm_1_25/', pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/perm_1_25/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_marine/perm_1_25/GLM_MARINe_chunk_*")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_perm_1_25", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(glm.model.output)
}

# ================================================================================== #

# Extract glm chunk for each environmental variable - perm 26:50

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/perm_26_50/', pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/perm_26_50/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_marine/perm_26_50/GLM_MARINe_chunk_*")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_perm_26_50", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(glm.model.output)
}

# ================================================================================== #

# Extract glm chunk for each environmental variable - perm 51:75

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/perm_51_75/', pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/perm_51_75/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_marine/perm_51_75/GLM_MARINe_chunk_*")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_perm_51_75", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(glm.model.output)
}

# ================================================================================== #

# Extract glm chunk for each environmental variable - perm 76:100

# Create list of file names for real data
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_marine/perm_76_100/', pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_marine/perm_76_100/', x))))

# Read all the files and add a column with the chunk
foreach(w=file_names_v, .errorhandling = "remove")%do%{  
    # State which file loading
    message(w)
    # Load file
    o = get(load(w))

    # Add column with identifier
    o %>% mutate(chunk = w) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_marine/perm_76_100/GLM_MARINe_chunk_*")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Chunk name
    c <- unique(tmp$chunk)

    # Extract data for a specific environmental variable
    tmp %>% filter(variable == env_var) -> tmp_env

    # Save subset
    save(tmp_env, file = paste(out_dir, "/glm_perm_76_100", env_var, "_", c, ".Rdata", sep = "") )

    # Clear load
    rm(glm.model.output)
}