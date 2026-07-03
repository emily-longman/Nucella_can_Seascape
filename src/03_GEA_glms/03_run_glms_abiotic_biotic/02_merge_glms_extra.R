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

# ================================================================================== #
# ================================================================================== #

# ================================================================================== #

# Real data
#load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")
# Perm data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.perm.Rdata")
glm.perm <- perm
# Load missing
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/GLM_chunk_missing.Rdata")

#real.missing <- glm.model.output[which(glm.model.output$data == "real"),]
#glm.real <- rbind(glm.real, real.missing)
#glm.real <- glm.real[order(match(glm.real$SNP_id, pooldata.snp.info$SNP_id)),]
#save(glm.real, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata"))


perm.missing <- glm.model.output[which(glm.model.output$data == "permutation"),]
glm.perm <- rbind(glm.perm, perm.missing)
glm.perm <- glm.perm[order(match(glm.perm$SNP_id, pooldata.snp.info$SNP_id)),]
save(glm.perm, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.perm.Rdata"))
