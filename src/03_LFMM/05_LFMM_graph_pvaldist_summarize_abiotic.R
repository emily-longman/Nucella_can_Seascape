# Graph and summarize abiotic LFMM output

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

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/LFMM_analyses")
if (!dir.exists(out_dir)) {dir.create(out_dir)}
out_fig_dir <- paste("output/figures/LFMM_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load abiotic real
load("data/processed/LFMM_analyses/LFMM_real_abiotic.RData")
str(LFMM_real)

# Extract real for abiotic var of interest
LFMM_real_var <- LFMM_real[which(LFMM_real$var == var),]

# Load perm
load(paste0("data/processed/LFMM_analyses/LFMM_perm_", var, ".RData"))
str(LFMM_perm)

# Join
LFMM <- rbind(LFMM_real_var, LFMM_perm)

# ================================================================================== #

# Bin p-vals
hist.obj.env = foreach(i = 0:max(LFMM$perm), .combine = "rbind")%do%{
    hist(LFMM$o.pvalues[LFMM$perm == i], breaks = 100) -> hist.obj
    data.frame(
    hist.obj$mids,
    hist.obj$counts,
    perm = i) -> o
}

# Graph pval distribution
pdf(paste0("output/figures/LFMM_analyses/LFMM_pval_dist_log_scale_", var, ".pdf"), width = 6.25, height = 5)
ggplot(hist.obj.env, aes(x=(hist.obj.mids), y=hist.obj.counts/1000, group=perm, color=perm==0, alpha=perm==0)) +
  geom_line(linewidth = 2) + 
  scale_alpha_manual(values = c("FALSE" = 0.1, "TRUE" = 1)) +
  scale_size_manual(values = c("FALSE" = 0.7, "TRUE" = 1.5)) +
  scale_color_manual(values = c("FALSE" = "#f63b3b",  "TRUE" = "black")) +
  scale_y_continuous(limits = c(0, 1500000), labels = function(x) format(x, scientific = FALSE)) +
  labs(title = paste0(var), x = "P-values", y = "# Thousand SNPs") +
  theme_bw(base_size = 30) +  theme(legend.position = "none") +
  scale_x_log10()
dev.off()

# ================================================================================== #

# Calc summary stats
#LFMM_summary <- foreach(i = 0:max(LFMM$perm), .combine = "rbind")%do%{
#    # Extract for just one perm
#    tmp <- LFMM %>% filter(perm == i)

#    # Create summary df
#    data.frame(
#        var = unique(tmp$var),
#        perm = i,
#        num_SNPs = dim(tmp)[1],
#        mean_zscores = mean(tmp$o.zscores),
#        num_SNPs_pval0.05 = sum(tmp$o.pvalues < 0.05),
#        mean_zscores_pval0.05 = mean(tmp$o.zscores[which(tmp$o.pvalues < 0.05)]),
#        num_SNPs_pval0.01 = sum(tmp$o.pvalues < 0.01),
#        mean_zscores_pval0.01 = mean(tmp$o.zscores[which(tmp$o.pvalues < 0.01)]),
#        num_SNPs_pval0.001 = sum(tmp$o.pvalues < 0.001),
#        mean_zscores_pval0.001 = mean(tmp$o.zscores[which(tmp$o.pvalues < 0.001)]),
#        num_SNPs_pvalbonferroni = sum(tmp$o.pvalues < 0.01/dim(tmp)[1]),
#        mean_zscores_pvalbonferroni = mean(tmp$o.zscores[which(tmp$o.pvalues < (0.01/dim(tmp)[1]))])
#    )
#}

# Save output
#save(LFMM_summary, file = paste0("data/processed/LFMM_analyses/LFMM_summary_", var, ".RData"))