# Merge glms output

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
#install.packages(c('rprojroot'))
#library(rprojroot)
# Specify root path
#root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
#setwd(root_path)

# ================================================================================== #

# Load packages
#install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach', 'pooldata'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)
library(pooldata)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/glms/cohensF")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Extract SNP info for all SNPs
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  pooldata.snp.info

# Rename columns
names(pooldata.snp.info)[1:2] = c("chr","pos")

# Make snp_id column
pooldata.snp.info <- pooldata.snp.info %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = "data/processed/glms/glms_chunk_analysis_abiotic/cohensF", pattern = "GLM_MARINe_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/glms/glms_chunk_analysis_abiotic/cohensF/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
glm.model.collated =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))

        # Add column with identifier
    o %>% mutate(chunk = i) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/glms/glms_chunk_analysis_abiotic/cohensF/GLM_MARINe_chunk_")) -> tmp

    # Remove end of chunk name
    tmp <- tmp %>% mutate(chunk = str_remove(chunk, pattern = ".Rdata"))

    # Return 
    tmp
}

# Check structure
str(glm.model.collated)

# ================================================================================== #

# Make snp_id column
glm.model.collated <- glm.model.collated %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

# Filter GLM to only sites in pooldata snp.info
glm.model.collated.filt <- glm.model.collated %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)

# Check structure
str(glm.model.collated.filt)

# ================================================================================== #

# Save merged data
save(glm.model.collated, file = "data/processed/glms/cohensF/glm_cohensF_filt_abiotic.Rdata")
