# Graph window analysis

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Load Data

# Read in outlier window
win.ph.mean.outliers <- read.csv("data/processed/GEA/glms/glms_window_summary/win.ph.mean.outliers.csv", header=T)

# Load GLM real data (glm.model.collated.filt.real)
load("data/processed/GEA/glms/glms_per_var/glm.collated_ph_mean_filt_real.Rdata")

# ================================================================================== #

# Identify SNPs in outlier windows
glm.ph.mean.win.outliers <- foreach(win.i=unique(win.ph.mean.outliers$win), .combine="rbind", .errorhandling="remove")%do%{
    
    # Extract window
    win.tmp <- win.ph.mean.outliers[which(win.ph.mean.outliers$win==win.i),]

    # Extract SNPs in window
    glm.tmp <- glm.model.collated.filt.real %>% filter(
        glm.model.collated.filt.real$chr == win.tmp$chr.x & 
        glm.model.collated.filt.real$pos > win.tmp$pos_min &
        glm.model.collated.filt.real$pos < win.tmp$pos_max)
}

# Remove duplicates
glm.ph.mean.win.outliers <- glm.ph.mean.win.outliers %>% distinct()

# ================================================================================== #

# Save
save(glm.ph.mean.win.outliers, file="data/processed/GEA/glms/glms_window_summary/glm.ph.mean.win.outliers.Rdata")