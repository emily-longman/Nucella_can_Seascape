# Use AIC to calc p-val

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
w = as.numeric(args[1]) # Chunk

# Prevent scientific notation
options(scipen=999)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load data

# Real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")
# Perm data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.perm.Rdata")

# Load snp info with groups
load("data/processed/GEA/glms/snp.guide.file.Rdata")

# ================================================================================== #

# Subset real and perm data based on chunk

# Subset snp.guide. file
snp.chunk <- snp.guide.file %>% filter(.groups == w)
# Subset real data
glm.real.chunk <- glm.real %>% filter(SNP_id %in% snp.chunk$SNP_id)
# Subset perm data
#glm.perm.chunk <- glm.perm %>% filter(SNP_id %in% snp.chunk$SNP_id)

# ================================================================================== #

# Identify which model is best for the real data

# Identify which column is the minimum
glm.real.chunk$minAIC <- names(glm.real.chunk[, 8:13])[apply(glm.real.chunk[, 8:13], MARGIN = 1, FUN = which.min)]
glm.real.chunk$minAIC_value <- do.call(pmin, c(glm.real.chunk[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.real.chunk <- glm.real.chunk %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# Check structure of real
str(glm.real.chunk)

# ================================================================================== #

# Identify which model is best for the perm data

# Identify which column is the minimum
glm.perm.chunk$minAIC_value <- do.call(pmin, c(glm.perm.chunk[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.perm.chunk <- glm.perm.chunk %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# ================================================================================== #

# Use permutations to calc p-val for SNPs where model with interaction is best fit 
o = foreach(i=1:dim(glm.real.chunk)[1], .combine = "rbind")%do%{
    
    # Extract real data for focal SNP
    real.tmp <- glm.real.chunk[i,]
    
    # Extract perm data for focal SNP
    perm.tmp <- glm.perm.chunk[which(glm.perm.chunk$SNP_id == real.tmp$SNP_id),]

    # Identify which column is the minimum
    #perm.tmp$minAIC_value <- do.call(pmin, c(perm.tmp[, 8:13], na.rm = TRUE))

    # Calc delta AIC between AIC_dem_both_int and model with min AIC
    #perm.tmp <- perm.tmp %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

    # Make data table and calculate p-val based on permutations (proportion of times that the perm data is less than the real data)
    data.frame(
          chr = unique(real.tmp$chr),
          pos = unique(real.tmp$pos),
          SNP_id = unique(real.tmp$SNP_id),
          real_minAIC = c(real.tmp$minAIC),
          real_minAIC_value = c(real.tmp$minAIC_value),
          real_deltaAIC = c(real.tmp$deltaAIC),
          p_val = c(mean(perm.tmp$deltaAIC <= real.tmp$deltaAIC)))
}

# Check structure
str(o)

# ================================================================================== #

# Generate folders and save output

# Folder name
folder_name <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/output_chunks")
if (!dir.exists(folder_name)) {dir.create(folder_name)}

# Save file for chunk w
file_name <- paste("GLM_output_chunk_", w, sep = "")
save(o, file = paste(folder_name, "/", file_name, ".Rdata", sep = ""))

# ================================================================================== #

message("done")
