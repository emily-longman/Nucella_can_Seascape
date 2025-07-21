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

out_dir <- paste("data/processed/GEA/glms/glms_output_by_model")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load
#load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")
load("data/processed/GEA/glms/glms_chunk_analysis/GLM_100perm_Bio-Oracle_chunk_1.Rdata")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Get names of enviro variables
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# ================================================================================== #

# Create SNP_id column
glm.model.output <- glm.model.output %>% mutate(SNP_id = paste(chr, pos, sep = "_"))





mod <- glm.model.output %>% group_by(variable)
