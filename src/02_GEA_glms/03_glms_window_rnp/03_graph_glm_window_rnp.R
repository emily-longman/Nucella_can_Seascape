# Create windows 

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

# Generate output directories

# Data directory
out_dir <- paste("data/processed/GEA/glms/glms_window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load GLM window rnp data
load("data/processed/GEA/glms/glms_window_summary/glm_windows_thetao_mean.RData")

# Load windows
load("data/processed/GEA/glms/glms_window_summary/windows.RData")

# ================================================================================== #

# Graph rnp p

# Create unique Chromosome number
chr.unique <- unique(wins_sum$chr)
wins_sum$chr.unique <- as.numeric(factor(wins_sum$chr, levels = chr.unique))

# Graph rnp
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_rnp_0.01_geompoint.pdf", width = 8, height = 8)
ggplot(wins_sum, aes(y=-log10(rnp.binom.p.0.01), x=chr.unique, col=factor(perm))) + 
  geom_point(alpha=0.8, size=1.3) + 
  scale_color_manual(values = c("black", rep("red", 100))) +
  facet_wrap(~variable) + ylim(0,1000) +
  theme_bw() + theme(legend.position = "none")
dev.off()

pdf("output/figures/GEA/glms/glms_window_summary/glm_window_rnp_0.01_geomline.pdf", width = 8, height = 8)
ggplot(wins_sum, aes(y=-log10(rnp.binom.p.0.01), x=chr.unique, col=factor(perm))) + 
  geom_line( ) + 
  scale_color_manual(values = c("black", rep("red", 100))) +
  facet_wrap(~variable) + ylim(0,1000) +
  theme_bw() + theme(legend.position = "none")
dev.off()

# ================================================================================== #
