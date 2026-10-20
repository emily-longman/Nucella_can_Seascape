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
out_fig_dir <- paste("output/figures/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Load SNP metadata

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# ================================================================================== #

# Read in Baypass  files

# Load BF output for 5 replicate Baypass runs
bf <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste0("data/processed/baypass/biotic/Mcali_", var, "/NC_biotic_Mcali_", var, "_run", i, "_summary_betai_reg.out"))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Change column names
setnames(bf, "BF(dB)", "bf_db")

# ================================================================================== #

# Average BF across replicate runs
bf.sum <- bf %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(bf_db.mean = mean(bf_db), bf_db.median = median(bf_db),
            bf_db.var=var(bf_db), eBPis.mean=mean(eBPis), eBPis.median=median(eBPis), eBPis.var=var(eBPis))


# Save
save(bf.sum, file=paste0("data/processed/baypass/biotic/bf.", sub("_", ".", var),".sum.Rdata"))
