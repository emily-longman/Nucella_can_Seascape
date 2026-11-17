# Merge annotations

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
#install.packages(c('data.table', 'tidyverse', 'dplyr', 'foreach', 'poolfstat'))
library(data.table)
library(tidyverse)
library(dplyr)
library(foreach)
library(poolfstat)

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/outlier_analyses/annotate_all")
file_names = as.list(dir(path = path, pattern = "chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/outlier_analyses/annotate_all/", sep = ""), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
o.merge =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# ================================================================================== #

# Save dataset - 14,897,468 SNPs
#save(o.merge, file="data/processed/outlier_analyses/Ncan.gds.annotations.RData")

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

# Filter for only SNPs in pooldata object
o.merge.filt <- o.merge %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)

# ================================================================================== #

# Save dataset - 8,191,999 SNPs
save(o.merge.filt, file="data/processed/outlier_analyses/Ncan.pooldata.annotations.RData")
