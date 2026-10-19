# Extract and merge LFMM output

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

#BiocManager::install("LEA")
library(LEA)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load 
load("data/processed/LFMM_analyses/LFMM_real_abiotic.RData")
str(LFMM_real)


# Bin data
hist.obj.env = foreach(i = unique(LFMM_real$var), .combine = "rbind")%do%{
hist(LFMM_real$o.pvalues, breaks = 100) -> hist.obj
data.frame(
  hist.obj$mids,
  hist.obj$counts,
  var = i) -> o
}

str(hist.obj.env)

# Graph pval distribution
pdf("output/figures/LFMM_analyses/LFMM_pval_dist_log_scale.pdf", width = 12, height = 8)
ggplot(hist.obj.env, aes(x=(hist.obj.mids), y=hist.obj.counts)) +
  geom_line(linewidth = 2) + 
  facet_wrap(~var, ncol = 4) +
  #scale_alpha_manual(values = c(0.1, 1)) +
  #scale_size_manual(values = c(0.7, 1.5)) +
  #scale_y_continuous(labels = function(x) format(x, scientific = FALSE)) +
  labs(x = "P-values", y = "Number of SNPs") +
  theme_bw(base_size = 30) +  theme(legend.position = "none") +
  scale_x_log10()
dev.off()