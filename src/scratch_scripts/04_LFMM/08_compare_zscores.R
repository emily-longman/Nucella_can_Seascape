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

# Mcali
load("data/processed/LFMM_analyses/LFMM_real_biotic_Mcali.RData")

# ================================================================================== #

# Calc summary stats
abiotic_summary <- foreach(i = unique(LFMM_real$var), .combine = "rbind")%do%{
    # Extract for just one perm
    tmp <- LFMM_real %>% filter(var == i)

    # Filter for NA
    tmp_filt <- tmp[!is.nan(tmp$o.zscores),]
    
    # Create summary df
    data.frame(
        var = unique(tmp_filt$var),
        num_SNPs_all = dim(tmp)[1],
        num_SNPs_filt = dim(tmp_filt)[1],
        num_SNPs_pval0.05 = sum(tmp_filt$o.pvalues < 0.05),
        mean_zscores_pval0.05 = mean(abs(tmp_filt$o.zscores[which(tmp_filt$o.pvalues < 0.05)])),
        num_SNPs_pval0.01 = sum(tmp_filt$o.pvalues < 0.01),
        mean_zscores_pval0.01 = mean(abs(tmp_filt$o.zscores[which(tmp_filt$o.pvalues < 0.01)])),
        num_SNPs_pval0.001 = sum(tmp_filt$o.pvalues < 0.001),
        mean_zscores_pval0.001 = mean(abs(tmp_filt$o.zscores[which(tmp_filt$o.pvalues < 0.001)])),
        num_SNPs_pvalbonferroni = sum(tmp_filt$o.pvalues < 0.01/dim(tmp_filt)[1]),
        mean_zscores_pvalbonferroni = mean(abs(tmp_filt$o.zscores[which(tmp_filt$o.pvalues < (0.01/dim(tmp_filt)[1]))]))
    )
}




