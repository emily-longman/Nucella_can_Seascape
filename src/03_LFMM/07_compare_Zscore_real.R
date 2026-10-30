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

# Load real data

# Abiotic
load("data/processed/LFMM_analyses/LFMM_real_abiotic.RData")
LFMM_real_abiotic <- LFMM_real

# Mcali
load("data/processed/LFMM_analyses/LFMM_real_biotic_Mcali.RData")
LFMM_real_Mcali <- LFMM_real

# MARINe
load("data/processed/LFMM_analyses/LFMM_real_biotic_MARINe.RData")
LFMM_real_MARINe <- LFMM_real

# ================================================================================== #

# Create hist object
breaks <- seq(0, 1, length.out = 101)
# Bin midpoints
midpoints <- head(breaks, -1) + diff(breaks) / 2

# ================================================================================== #

# List of variables
vars_abiotic <- unique(LFMM_real_abiotic$var)

# Calc summary stats
summary_abiotic <- foreach(i = vars_abiotic, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_abiotic %>% filter(var == i)
    
    # Bin pvalues
    tmp$bin <- cut(tmp$o.pvalues, breaks = breaks, include.lowest = TRUE)
   
    # Summarize
    o <- tmp %>% group_by(bin) %>%
        summarize(counts = n(), mean_zscore = mean(o.zscores, na.rm = TRUE), mean_abs_zscore = mean(abs(o.zscores), na.rm = TRUE))
    
    # Add midpoint
    o$mids <- midpoints
    o$var <- i

    # Return
    o
}

# Graph zscore vs p-values
pdf("output/figures/LFMM_analyses/LFMM_pval_zscores_abiotic.pdf", width = 14, height = 8)
ggplot(summary_abiotic, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    theme_bw(base_size = 30)
dev.off()

# Graph zscores vs p-values (log scale)
pdf("output/figures/LFMM_analyses/LFMM_log_pval_zscores_abiotic.pdf", width = 14, height = 8)
ggplot(summary_abiotic, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    scale_x_log10() +
    theme_bw(base_size = 30)
dev.off()

###########

# Calc summary stats
summary_abiotic_sigSNPs <- foreach(i = vars_abiotic, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_abiotic %>% filter(var == i)
    # Summarize
    q <- quantile(abs(tmp$o.zscores), probs = c(0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i, 
        n_sig = sum(tmp$o.pvalues < 0.05),
        mean_abs_zscore_0.05 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        median_abs_zscore_0.05 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        mean_abs_zscore_0.01 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        median_abs_zscore_0.01 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        mean_abs_zscore_0.001 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        median_abs_zscore_0.001 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        max_abs_zscore =  max(abs(tmp$o.zscores)),
        abs_zscore_quant_95 = q[1],
        abs_zscore_quant_99 = q[2],
        abs_zscore_quant_99.5 = q[3],
        abs_zscore_quant_99.9 = q[4]
    )
}

# Join with full names
summary_abiotic_sigSNPs <- left_join(summary_abiotic_sigSNPs, var_full_names, by="variable")

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_abiotic.pdf", width = 9, height = 5.5)
ggplot(summary_abiotic_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for p<0.05
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.05_abiotic.pdf", width = 9, height = 5.5)
ggplot(summary_abiotic_sigSNPs, aes(x = mean_abs_zscore_0.05, y = reorder(variable_full_name_edit, mean_abs_zscore_0.05))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.05"), y = "") +
  xlim(0, 5.5) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_abiotic.pdf", width = 9, height = 5.5)
ggplot(summary_abiotic_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  xlim(0, 10.5) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #

# List of variables
vars_Mcali <- unique(LFMM_real_Mcali$var)

# Calc summary stats
summary_biotic_Mcali <- foreach(i = vars_Mcali, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_Mcali %>% filter(var == i)
    
    # Bin pvalues
    tmp$bin <- cut(tmp$o.pvalues, breaks = breaks, include.lowest = TRUE)
   
    # Summarize
    o <- tmp %>% group_by(bin) %>%
        summarize(counts = n(), mean_zscore = mean(o.zscores, na.rm = TRUE), mean_abs_zscore = mean(abs(o.zscores), na.rm = TRUE))
    
    # Add midpoint
    o$mids <- midpoints
    o$var <- i

    # Return
    o
}

# Graph zscore vs p-values
pdf("output/figures/LFMM_analyses/LFMM_pval_zscores_biotic_Mcali.pdf", width = 14, height = 8)
ggplot(summary_biotic_Mcali, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    theme_bw(base_size = 30)
dev.off()

# Graph zscores vs p-values (log scale)
pdf("output/figures/LFMM_analyses/LFMM_log_pval_zscores_biotic_Mcali.pdf", width = 14, height = 8)
ggplot(summary_biotic_Mcali, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    scale_x_log10() +
    theme_bw(base_size = 30)
dev.off()


###########

# Calc summary stats
summary_biotic_Mcali_sigSNPs <- foreach(i = vars_Mcali, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_Mcali %>% filter(var == i)
    # Summarize
    q <- quantile(abs(tmp$o.zscores), probs = c(0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i, 
        n_sig = sum(tmp$o.pvalues < 0.05),
        mean_abs_zscore_0.05 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        median_abs_zscore_0.05 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        mean_abs_zscore_0.01 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        median_abs_zscore_0.01 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        mean_abs_zscore_0.001 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        median_abs_zscore_0.001 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        max_abs_zscore =  max(abs(tmp$o.zscores)),
        abs_zscore_quant_95 = q[1],
        abs_zscore_quant_99 = q[2],
        abs_zscore_quant_99.5 = q[3],
        abs_zscore_quant_99.9 = q[4]
    )
}

# Join with full names
summary_biotic_Mcali_sigSNPs <- left_join(summary_biotic_Mcali_sigSNPs, var_full_names, by="variable")

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_biotic_Mcali.pdf", width = 9, height = 5.5)
ggplot(summary_biotic_Mcali_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for p<0.05
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.05_biotic_Mcali.pdf", width = 9, height = 3.5)
ggplot(summary_biotic_Mcali_sigSNPs, aes(x = mean_abs_zscore_0.05, y = reorder(variable_full_name_edit, mean_abs_zscore_0.05))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.05"), y = "") +
  xlim(0, 4.5) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_biotic_Mcali.pdf", width = 9, height = 3.5)
ggplot(summary_biotic_Mcali_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  xlim(0, 7) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #

# List of variables
vars_MARINe <- unique(LFMM_real_MARINe$var)

# Calc summary stats
summary_biotic_MARINe <- foreach(i = vars_MARINe, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_MARINe %>% filter(var == i)
    
    # Bin pvalues
    tmp$bin <- cut(tmp$o.pvalues, breaks = breaks, include.lowest = TRUE)
   
    # Summarize
    o <- tmp %>% group_by(bin) %>%
        summarize(counts = n(), mean_zscore = mean(o.zscores, na.rm = TRUE), mean_abs_zscore = mean(abs(o.zscores), na.rm = TRUE))
    
    # Add midpoint
    o$mids <- midpoints
    o$var <- i

    # Return
    o
}

# Graph zscore vs p-values
pdf("output/figures/LFMM_analyses/LFMM_pval_zscores_biotic_MARINe.pdf", width = 14, height = 8)
ggplot(summary_biotic_MARINe, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    theme_bw(base_size = 30)
dev.off()

# Graph zscores vs p-values (log scale)
pdf("output/figures/LFMM_analyses/LFMM_log_pval_zscores_biotic_MARINe.pdf", width = 14, height = 8)
ggplot(summary_biotic_MARINe, aes(x=mids, y=mean_abs_zscore, color = var)) +
    geom_line() + 
    labs(x = "P-val bin", y = "Mean Z-score") +
    scale_x_log10() +
    theme_bw(base_size = 30)
dev.off()

###########

# Calc summary stats
summary_biotic_MARINe_sigSNPs <- foreach(i = vars_MARINe, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_real_MARINe %>% filter(var == i)
    # Summarize
    q <- quantile(abs(tmp$o.zscores), probs = c(0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i, 
        n_sig = sum(tmp$o.pvalues < 0.05),
        mean_abs_zscore_0.05 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        median_abs_zscore_0.05 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        mean_abs_zscore_0.01 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        median_abs_zscore_0.01 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        mean_abs_zscore_0.001 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        median_abs_zscore_0.001 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        max_abs_zscore =  max(abs(tmp$o.zscores)),
        abs_zscore_quant_95 = q[1],
        abs_zscore_quant_99 = q[2],
        abs_zscore_quant_99.5 = q[3],
        abs_zscore_quant_99.9 = q[4]
    )
}

# Join with full names
summary_biotic_MARINe_sigSNPs <- left_join(summary_biotic_MARINe_sigSNPs, var_full_names, by="variable")

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_biotic_MARINe.pdf", width = 9, height = 5.5)
ggplot(summary_biotic_MARINe_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for p<0.05
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.05_biotic_MARINe.pdf", width = 9, height = 6.5)
ggplot(summary_biotic_MARINe_sigSNPs, aes(x = mean_abs_zscore_0.05, y = reorder(variable_full_name_edit, mean_abs_zscore_0.05))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.05"), y = "") +
  xlim(0, 3) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_biotic_MARINe.pdf", width = 9, height = 6.5)
ggplot(summary_biotic_MARINe_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  xlim(0, 20) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #

# Join biotic datasets
summary_biotic_sigSNPs <- rbind(summary_biotic_Mcali_sigSNPs, summary_biotic_MARINe_sigSNPs)

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_biotic.pdf", width = 9, height = 9)
ggplot(summary_biotic_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for p<0.05
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.05_biotic.pdf", width = 9, height = 9)
ggplot(summary_biotic_sigSNPs, aes(x = mean_abs_zscore_0.05, y = reorder(variable_full_name_edit, mean_abs_zscore_0.05))) + 
  geom_point(size = 8) + 
  labs(x = expression("Mean " * abs(z)* " for p<0.05"), y = "") +
  xlim(0, 4.5) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_biotic.pdf", width = 9, height = 9)
ggplot(summary_biotic_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  xlim(0, 20) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #

# Join full dataset
LFMM_summary <- rbind(summary_abiotic_sigSNPs, summary_biotic_sigSNPs)
# Save
save(LFMM_summary, file = "data/processed/LFMM_analyses/LFMM_summary.RData")