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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'magrittr', 'ggtext'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(magrittr)
library(ggtext)

# ================================================================================== #

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)

# Load LFMM results
load("data/processed/LFMM_analyses/LFMM_all.RData")
str(LFMM_all)

# ================================================================================== #

# List of variables
vars <- unique(LFMM_all$var)

# Calc summary stats
summary_all_sigSNPs <- foreach(i = vars, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- LFMM_all %>% filter(var == i)
    # Summarize
    q <- quantile(abs(tmp$o.zscores), probs = c(0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i, 
        n_tot = dim(tmp)[1],
        n_sig_0.05 = sum(tmp$o.pvalues < 0.05),
        n_sig_0.01 = sum(tmp$o.pvalues < 0.01),
        n_sig_0.001 = sum(tmp$o.pvalues < 0.001),
        mean_abs_zscore_0.05 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        median_abs_zscore_0.05 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.05)])),
        mean_abs_zscore_0.01 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        median_abs_zscore_0.01 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.01)])),
        mean_abs_zscore_0.001 = mean(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        median_abs_zscore_0.001 = median(abs(tmp$o.zscores[which(tmp$o.pvalues < 0.001)])),
        abs_zscore_quant_95 = q[1],
        abs_zscore_quant_99 = q[2],
        abs_zscore_quant_99.5 = q[3],
        abs_zscore_quant_99.9 = q[4]
    )
}

# Join with full names
summary_sigSNPs <- left_join(summary_all_sigSNPs, var_full_names, by="variable")

# Write output
write.csv(summary_sigSNPs, "data/processed/LFMM_analyses/summary_sigSNPs.csv")

# ================================================================================== #

# For graphing make species names italics
summary_sigSNPs$variable_full_name_italics <- 
c("SST (max)", "SST (min)", "SST (range)", "SST (mean)", 
"Chlorophyll (mean)", "O2 (mean)", "pH (min)", "pH (mean)", "Salinity (mean)",
"*B. glandula* (mean)", "*B. glandula* (hmean)", "*C. dalli* (mean)", "*C. dalli* (hmean)",
"*M. californianus* (mean)", "*M. californianus* (hmean)", "*M. trossulus* (mean)", "*M. trossulus* (hmean)", 
"*N. canaliculata* (mean)", "*N. canaliculata* (hmean)","*N. ostrina* (mean)", "*N. ostrina* (hmean)", 
"*P. ochraceus* (mean)", "*P. ochraceus* (hmean)",  "STI (mean)", "Cross-sectional thk (mean)", "Max thk (mean)", "Min thk (mean)")

# Graph LFMM results

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_all.pdf", width = 12, height = 12)
ggplot(summary_sigSNPs, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_italics, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  xlim(0.4, 9) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_markdown(size = 22))
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_all.pdf", width = 12, height = 12)
ggplot(summary_sigSNPs, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_italics, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  xlim(0.85, 19.5) +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_markdown(size = 22))
dev.off()


# Split up by abiotic and biotic
summary_sigSNPs_abiotic <- summary_sigSNPs %>% filter(Variable == "Abiotic")
summary_sigSNPs_biotic <- summary_sigSNPs %>% filter(Variable == "Biotic")

# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_all_abiotic.pdf", width = 8, height = 5)
ggplot(summary_sigSNPs_abiotic, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22)) 
dev.off()
# Graph z scores for p<0.001
pdf("output/figures/LFMM_analyses/LFMM_zscore_pval0.001_all_biotic.pdf", width = 9, height = 10)
ggplot(summary_sigSNPs_biotic, aes(x = mean_abs_zscore_0.001, y = reorder(variable_full_name_edit, mean_abs_zscore_0.001))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean " * abs(z)* " for p<0.001"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22)) 
dev.off()

# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_all.pdf", width = 8, height = 5)
ggplot(summary_sigSNPs_abiotic, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()
# Graph z scores for 99.9%
pdf("output/figures/LFMM_analyses/LFMM_zscore_99.9perc_all.pdf", width = 9, height = 10)
ggplot(summary_sigSNPs_biotic, aes(x = abs_zscore_quant_99.9, y = reorder(variable_full_name_edit, abs_zscore_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile " * abs(z)), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()
