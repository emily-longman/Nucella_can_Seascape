# Compare outliers for GLM and Baypass

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

# Load GLM outliers (glm.ph.mean.win.outliers)
load("data/processed/GEA/glms/glms_window_summary/glm.ph.mean.win.outliers.Rdata")

# Load Baypass outliers
bf.ph.mean.sum.outliers <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.csv", header=T)

# ================================================================================== #

# Identify overlap - 173 SNPs
common_values <- intersect(glm.ph.mean.win.outliers$SNP_id, bf.ph.mean.sum.outliers$SNP_id)

# Extract from each
glm.overlap <- glm.ph.mean.win.outliers %>% filter(SNP_id %in% common_values)
baypass.overlap <- bf.ph.mean.sum.outliers %>% filter(SNP_id %in% common_values)

# Merge
ph.mean.outliers <- left_join(glm.overlap, baypass.overlap, by=c("chr", "pos", "SNP_id"))

# ================================================================================== #

# Save
save(ph.mean.outliers, file="data/processed/outlier_analyses/ph.mean.outliers.Rdata")
write.csv(ph.mean.outliers, "data/processed/outlier_analyses/ph.mean.outliers.csv", row.names = F, quote = F)