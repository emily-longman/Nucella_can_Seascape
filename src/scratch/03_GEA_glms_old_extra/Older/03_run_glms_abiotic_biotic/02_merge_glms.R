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

# Read all the files and add a column with the chunk
glm.real =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.real)

# Save merged data
save(glm.real, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata"))

# ================================================================================== #
# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/perm_1_50/', pattern = "GLM_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/perm_1_50/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
glm.perm1.50 =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.perm1.50)

# ================================================================================== #
# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/perm_51_100/', pattern = "GLM_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/perm_51_100/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
glm.perm51.100 =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.perm51.100)

# Save merged data
#save(glm.model.collated.perm51.100, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.perm51.100.Rdata"))

# ================================================================================== #
# ================================================================================== #

# Join the perm data
glm.perm <- rbind(glm.perm1.50, glm.perm51.100)

# Save merged data
save(glm.perm, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.perm.Rdata"))


# ================================================================================== #

# Real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")
# Perm data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.perm.Rdata")
glm.perm <- perm
# Load missing
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/GLM_chunk_missing.Rdata")

real.missing <- glm.model.output[which(glm.model.output$data == "real"),]
glm.real <- rbind(glm.real, real.missing)
glm.real <- glm.real[order(match(glm.real$SNP_id, pooldata.snp.info$SNP_id)),]
save(glm.real, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata"))


perm.missing <- glm.model.output[which(glm.model.output$data == "permutation"),]
glm.perm <- rbind(glm.perm, perm.missing)
glm.perm <- glm.perm[order(match(glm.perm$SNP_id, pooldata.snp.info$SNP_id)),]
save(glm.perm, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.perm.Rdata"))
