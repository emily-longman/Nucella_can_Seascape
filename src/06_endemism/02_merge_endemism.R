# Merge endemism analysis
# Modified from L. Proud; Original from A. Bergland

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
#install.packages(c('SeqArray', 'data.table', 'tidyverse', 'dplyr', 'geodist', 'foreach'))
library(SeqArray)
library(data.table)
library(tidyverse)
library(dplyr)
library(geodist)
library(foreach)

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/endemism/chunks")
file_names = as.list(dir(path = path, pattern = "slice_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/endemism/chunks/", sep = ""), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
endemism.merge =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# ================================================================================== #

# Save dataset - 14,526,366 SNPs (although GDS is 14,897,468 SNPs)
save(endemism.merge, file="data/processed/endemism/endemism.merge.RData")

# ================================================================================== #

# Load pooldata object - 8,277,206 SNPs
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

# Filter for only SNPs in pooldata object - 8,191,999 SNPs
endemism.merge.filt <- endemism.merge %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)