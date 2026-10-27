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

# Create list of file names
file_names = as.list(dir(path = 'data/processed/LFMM_analyses/', pattern = paste0("LFMM_summary_*")))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/LFMM_analyses/', x))))

# Read all the files
LFMM_sum <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Format
LFMM_sum$var <- as.factor(LFMM_sum$var)
# Filter NAs
LFMM_sum_filt <- LFMM_sum[!is.nan(LFMM_sum$mean_zscores_pvalbonferroni),]

# Graph bonferroni zscores
pdf("output/figures/LFMM_analyses/LFMM_filt_zscore_abiotic.pdf", width = 7, height = 8)
ggplot(LFMM_sum_filt, aes(x=mean_zscores_pvalbonferroni, y=var)) +
    geom_boxplot(data = LFMM_sum_filt[which(LFMM_sum_filt$perm != 0),], aes(x=mean_zscores_pvalbonferroni, y=var), color = "grey") +
    geom_point(data = LFMM_sum_filt[which(LFMM_sum_filt$perm == 0),], aes(x=mean_zscores_pvalbonferroni, y=var), color = "black") + 
    labs(x = "Mean Z-score\n(Bonferroni)", y = "") +
    theme_bw(base_size = 30)
dev.off()
