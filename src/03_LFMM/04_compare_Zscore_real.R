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

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)

# Load LFMM results

# ================================================================================== #

# List of variables
vars <- unique(var_full_names$variable)

# Calc summary stats
summary_all_sigSNPs <- foreach(i = vars, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_all %>% filter(var == i)
    # Summarize
    q <- quantile(abs(tmp$o.zscores.mean), probs = c(0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i, 
        n_tot = dim(tmp)[1],
        n_sig_0.05 = sum(tmp$o.pvalues.mean < 0.05),
        n_sig_0.01 = sum(tmp$o.pvalues.mean < 0.01),
        n_sig_0.001 = sum(tmp$o.pvalues.mean < 0.001),
        mean_abs_zscore_0.05 = mean(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.05)])),
        median_abs_zscore_0.05 = median(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.05)])),
        mean_abs_zscore_0.01 = mean(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.01)])),
        median_abs_zscore_0.01 = median(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.01)])),
        mean_abs_zscore_0.001 = mean(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.001)])),
        median_abs_zscore_0.001 = median(abs(tmp$o.zscores.mean[which(tmp$o.pvalues.mean < 0.001)])),
        max_abs_zscore =  max(abs(tmp$o.zscores.mean)),
        abs_zscore_quant_95 = q[1],
        abs_zscore_quant_99 = q[2],
        abs_zscore_quant_99.5 = q[3],
        abs_zscore_quant_99.9 = q[4]
    )
}

# Join with full names
summary_sigSNPs <- left_join(summary_all_sigSNPs, var_full_names, by="variable")

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_all_5run.pdf", width = 9, height = 12)
ggplot(summary_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22)) 
dev.off()

# Graph z scores for p<0.05
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.05_all_5run.pdf", width = 9, height = 12)
ggplot(summary_sigSNPs, aes(x = mean_abs_zscore_0.05, y = reorder(variable_full_name_edit, mean_abs_zscore_0.05))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean " * abs(z)* " for p<0.05"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_all_5run.pdf", width = 9, height = 12)
ggplot(summary_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()