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

# Generate output directories

# Figure directory
out_data_dir <- paste("data/processed/baypass/enrichment")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Read in Baypass files

# Load BF output for 5 replicate Baypass runs
bf_real <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste0("data/processed/baypass/biotic/Mcali_", var, "/NC_biotic_Mcali_", var, "_run", i, "_summary_betai_reg.out"))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Change column names
setnames(bf_real, "BF(dB)", "bf_db")

# Average BF across replicate runs
bf_real_sum <- bf_real %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(bf_db.mean = mean(bf_db), bf_db.median = median(bf_db),
            bf_db.var=var(bf_db), eBPis.mean=mean(eBPis), eBPis.median=median(eBPis), eBPis.var=var(eBPis))

# ================================================================================== #

# Read in Baypass POD files

# Load BF output for 5 replicate Baypass runs
bf_POD <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste0("data/processed/baypass/biotic/Mcali_", var, "_POD/NC_biotic_Mcali_", var, "_POD_run", i, "_summary_betai_reg.out"))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Change column names
setnames(bf_POD, "BF(dB)", "bf_db")

# Calculate quantiles for each POD
bf_5POD_quant <- bf_POD %>% group_by(rep) %>% reframe(bf_db = quantile(bf_db, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()

# Average quantiles across POD runs
bf_5POD_thr <- bf_5POD_quant %>% group_by(thr) %>% summarize(bf_db.mean=mean(bf_db))

# Save
save(bf_5POD_thr, file=paste0("data/processed/baypass/biotic/bf_", var, "_5POD_thr.Rdata"))

# ================================================================================== #

# Enrichment

# Fold enrichment above POD threshold
bf_enrich <- bf_real_sum %>%
    summarize(
    var = var,
    n = dim(bf_real_sum)[1],
    count_95POD = sum(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.95)]),
    count_99POD = sum(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.99)]),
    count_99.9POD = sum(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.999)]),
    mean_BF_95POD_sig = mean(bf_db.mean[which(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.95)])]),
    mean_BF_99POD_sig = mean(bf_db.mean[which(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.99)])]),
    mean_BF_99.9POD_sig = mean(bf_db.mean[which(bf_db.mean > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.999)])]),
    BF_quant_95 = quantile(bf_db.mean, 0.95),
    BF_quant_99 = quantile(bf_db.mean, 0.99),
    BF_quant_99.9 = quantile(bf_db.mean, 0.999)
    ) %>% as.data.frame()

# Save
save(bf_enrich, file=paste0("data/processed/baypass/enrichment/bf_", var, "_enrichment.Rdata"))

# ================================================================================== E#

# Fold enrichment
bf_POD_enrichment <- bf_POD %>% group_by(rep) %>% 
    summarize(BF_10 = sum(bf_db > 10), BF_15 = sum(bf_db > 15), BF_20 = sum(bf_db > 20), BF_95POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.95)]), BF_99POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.99)]), BF_99.9POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.999)])) %>%
    summarize(BF_10 = mean(BF_10), BF_15 = mean(BF_15), BF_20 = mean(BF_20), BF_95POD = mean(BF_95POD), BF_99POD = mean(BF_99POD), BF_99.9POD = mean(BF_99.9POD))

# Real data per rep
bf_real_enrichment <- bf_real %>% group_by(rep) %>%
    summarize(
        var = var,
        E_BF_10 = sum(bf_db > 10)/bf_POD_enrichment$BF_10,
        E_BF_15 = sum(bf_db > 15)/bf_POD_enrichment$BF_15,
        E_BF_20 = sum(bf_db > 20)/bf_POD_enrichment$BF_20,
        E_BF_95POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.95)])/bf_POD_enrichment$BF_95POD,
        E_BF_99POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.99)])/bf_POD_enrichment$BF_99POD,
        E_BF_999POD = sum(bf_db > bf_5POD_thr$bf_db.mean[which(bf_5POD_thr$thr == 0.999)])/bf_POD_enrichment$BF_99.9POD,   
        )

# Save
save(bf_real_enrichment, file=paste0("data/processed/baypass/enrichment/bf_", var,"_fold_enrichment.Rdata"))
