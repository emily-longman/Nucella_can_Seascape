# Compare Z scores

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'magrittr'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(magrittr)

# ================================================================================== #

# Load real data

# Abiotic
LFMM_real <- load("data/processed/LFMM_analyses/LFMM_real_abiotic.RData")
LFMM_real_abiotic <- LFMM_real

# MARINe
LFMM_real <- load("data/processed/LFMM_analyses/LFMM_real_biotic_MARINe.RData")
LFMM_real_MARINe <- LFMM_real

# Mcali 
LFMM_real <- load("data/processed/LFMM_analyses/LFMM_real_biotic_Mcali.RData")
LFMM_real_Mcali <- LFMM_real

# ================================================================================== #

# Read all the files

# Graph bonferroni zscores
pdf("output/figures/LFMM_analyses/LFMM_filt_zscore_abiotic.pdf", width = 7, height = 8)
ggplot(LFMM_sum_filt, aes(x=mean_zscores_pvalbonferroni, y=var)) +
    geom_boxplot(data = LFMM_sum_filt[which(LFMM_sum_filt$perm != 0),], aes(x=mean_zscores_pvalbonferroni, y=var), color = "grey") +
    geom_point(data = LFMM_sum_filt[which(LFMM_sum_filt$perm == 0),], aes(x=mean_zscores_pvalbonferroni, y=var), color = "black") + 
    labs(x = "Mean Z-score\n(Bonferroni)", y = "") +
    theme_bw(base_size = 30)
dev.off()
