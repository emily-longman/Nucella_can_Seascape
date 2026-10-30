# Summarize Baypass

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

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)

# ================================================================================== #

# Read in  files

# File names
file_names = as.list(dir(path = 'data/processed/baypass/enrichment/', pattern = "*_enrichment.Rdata"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/enrichment/', x))))

# Load BF output for 5 replicate Baypass runs
bf_enrich_all <- foreach(i=file_names_v[seq(from = 1, to = length(file_names_v), by = 2) ], .combine = rbind)%do%{
    message(i)
    load(i)
    return(bf_enrich)
}

# Join
bf_enrich_all <- left_join(bf_enrich_all, var_full_names, by = "var")
bf_enrich_all$Variable <- as.factor(bf_enrich_all$Variable)

# ================================================================================== #

# Graph
pdf("output/figures/baypass/Baypass_BF_count_95POD.pdf", width = 9, height = 12)
ggplot(bf_enrich_all, aes(x = count_95POD, y = reorder(variable_full_name_edit, count_95POD))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 8) + 
    labs(x = "# SNPs > 95%POD Thr", y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_BF_count_99.9POD.pdf", width = 9, height = 12)
ggplot(bf_enrich_all, aes(x = count_99.9POD, y = reorder(variable_full_name_edit, count_99.9POD))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 8) + 
    labs(x = "# SNPs > 99.9%POD Thr", y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_mean_BF_99.9POD_sig.pdf", width = 10, height = 12)
ggplot(bf_enrich_all, aes(x = mean_BF_99.9POD_sig, y = reorder(variable_full_name_edit, mean_BF_99.9POD_sig))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 8) + 
    labs(x = "mean(BF) for POD sig SNPs", y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_BF_quant_99.9.pdf", width = 9, height = 12)
ggplot(bf_enrich_all, aes(x = BF_quant_99.9, y = reorder(variable_full_name_edit, BF_quant_99.9))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 8) + 
    labs(x = expression("BF 99.9% Percentile"), y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Compare baypass with LFMM
BayPass_LFMM_summary <- left_join(bf_enrich_all, LFMM_summary)
BayPass_LFMM_summary$Variable <- as.factor(BayPass_LFMM_summary$Variable)
BayPass_LFMM_summary$shape <- as.factor(c(rep("C", 5), "T", rep("C", 12), "T",rep("C", 7)))

# Graph
pdf("output/figures/Baypass_LFMM_quant_99.9.pdf", width = 7, height = 6)
ggplot(BayPass_LFMM_summary, aes(x = BF_quant_99.9, y = abs_zscore_quant_99.9, color = Variable, shape = shape)) + 
    geom_point(size = 8) + 
    scale_color_manual(values = c("#d2931f", "purple")) + 
    labs(x = expression("BF 99.9% Percentile"), y = expression("99.9% Percentile " * abs(z))) + 
    guides(shape = "none") +
    scale_x_continuous(limits = c(5,14)) + scale_y_continuous(limits = c(0,20)) +
    theme_linedraw(base_size=30) + theme( legend.position = c(0.8, 0.75), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5,), legend.title = element_blank())
dev.off()
pdf("output/figures/Baypass_LFMM_quant_99.9.pdf", width = 7, height = 6)
ggplot(BayPass_LFMM_summary, aes(y = BF_quant_99.9, x = abs_zscore_quant_99.9, color = Variable, shape = shape)) + 
    geom_point(size = 8) + 
    scale_color_manual(values = c("#d2931f", "purple")) + 
    labs(y = expression("BF 99.9% Percentile"), x = expression("99.9% Percentile " * abs(z))) + 
    guides(shape = "none") +
    scale_y_continuous(limits = c(5,14)) + scale_x_continuous(limits = c(0,20)) +
    theme_linedraw(base_size=30) + theme( legend.position = c(0.8, 0.75), legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5,), legend.title = element_blank())
dev.off()


# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Read in fold enrichment files

# Load BF output for 5 replicate Baypass runs
bf_fold_enrich_all <- foreach(i=file_names_v[seq(from = 2, to = length(file_names_v), by = 2) ], .combine = rbind)%do%{
    message(i)
    load(i)
    return(bf_real_enrichment)
}

# Join
bf_fold_enrich_all <- left_join(bf_fold_enrich_all, var_full_names, by = "var")
bf_fold_enrich_all$Variable <- as.factor(bf_fold_enrich_all$Variable)

# Graph
pdf("output/figures/baypass/Baypass_BF_fold_enrich_10BF.pdf", width = 9, height = 12)
ggplot(bf_fold_enrich_all, aes(x = E_BF_10, y = reorder(variable_full_name_edit, E_BF_10))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 6, alpha = 0.8) + 
    labs(x = expression("Fold Enrich for BF>10"), y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_BF_fold_enrich_15BF.pdf", width = 9, height = 12)
ggplot(bf_fold_enrich_all, aes(x = E_BF_15, y = reorder(variable_full_name_edit, E_BF_15))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 6, alpha = 0.8) + 
    labs(x = expression("Fold Enrich for BF>15"), y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_BF_fold_enrich_20BF.pdf", width = 9, height = 12)
ggplot(bf_fold_enrich_all, aes(x = E_BF_20, y = reorder(variable_full_name_edit, E_BF_20))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 6, alpha = 0.8) + 
    labs(x = expression("Fold Enrich for BF>20"), y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph
pdf("output/figures/baypass/Baypass_BF_fold_enrich_BF_999POD.pdf", width = 9, height = 12)
ggplot(bf_fold_enrich_all, aes(x = E_BF_999POD, y = reorder(variable_full_name_edit, E_BF_999POD))) + 
    facet_wrap(~Variable, ncol = 1, scales = "free_y") + 
    geom_point(size = 6, alpha = 0.8) + 
    labs(x = expression("Fold Enrich 99.9% POD"), y = "") +
    theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()
