# Create guide file

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'groupdata2'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(groupdata2)

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

# Group the snp list into 920 chunks with 9,000 SNPs in each
snp.guide.file <- group(pooldata.snp.info, n=9000, method = "greedy")

# ================================================================================== #

# Save
save(snp.guide.file, file = paste("data/processed/GEA/glms/snp.guide.file.Rdata"))