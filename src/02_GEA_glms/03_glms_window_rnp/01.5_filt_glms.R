# Filt glms output

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
#install.packages(c('tidyverse', 'dplyr'))
library(tidyverse)
library(dplyr)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
env_var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# State variable name
message(env_var)

# Load data
load(paste0("data/processed/GEA/glms/glms_per_env_var/glm.collated_", env_var, ".Rdata") )

# Load data
#load("data/processed/GEA/glms/glms_per_env_var/glm.collated_ph_mean.Rdata")

# Make snp_id column
glm.model.collated <- glm.model.collated %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

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

# Filter GLM to only sites in pooldata snp.info
glm.model.collated.filt <- glm.model.collated %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)

# Check structure
str(glm.model.collated.filt)

# ================================================================================== #

# Save glm output
#save(glm.model.collated.filt, file="data/processed/GEA/glms/glms_per_env_var/glm.collated_ph_mean_8M.Rdata")
save(glm.model.collated.filt, file=paste("data/processed/GEA/glms/glms_per_env_var/glm.collated_", env_var, "_8M.Rdata", sep = ""))
