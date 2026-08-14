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

# Load GLM outliers (glm.Mtross.mean.win.outliers)
load("data/processed/GEA/glms/glms_window_summary/glm.Mtross.mean.win.outliers.Rdata")

# Load Baypass outliers
bf.Mtross.mean.sum.outliers <- read.csv("data/processed/baypass/bf.Mtross.mean.sum.outliers.csv", header=T)

# ================================================================================== #

# Identify overlap - 173 SNPs
common_values <- intersect(glm.Mtross.mean.win.outliers$SNP_id, bf.Mtross.mean.sum.outliers$SNP_id)

# Extract from each
glm.overlap <- glm.Mtross.mean.win.outliers %>% filter(SNP_id %in% common_values)
baypass.overlap <- bf.Mtross.mean.sum.outliers %>% filter(SNP_id %in% common_values)

# Merge
Mtross.mean.outliers <- left_join(glm.overlap, baypass.overlap, by=c("chr", "pos", "SNP_id"))

# ================================================================================== #

# Save
save(Mtross.mean.outliers, file="data/processed/outlier_analyses/Mtross.mean.outliers.Rdata")
write.csv(Mtross.mean.outliers, "data/processed/outlier_analyses/Mtross.mean.outliers.csv", row.names = F, quote = F)