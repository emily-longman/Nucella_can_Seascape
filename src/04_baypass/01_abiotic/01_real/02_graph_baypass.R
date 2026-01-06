# Graph window analysis

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
out_fig_dir <- paste("output/figures/GEA/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load Data

# Read in SNP data
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")


# Read in Baypass results
baypass_betai <- read.table("data/processed/GEA/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_summary_betai.out", header=T)

# Join baypass results with snp metadata
baypass_betai_pos <- cbind(snp.meta, baypass_betai)


# ================================================================================== #

# Graph rnp p

# Create unique Chromosome number
chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = chr.unique))

# Graph rnp geompoint
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_rnp_0.01_geompoint.pdf", width = 8, height = 8)
plot(baypass_betai_pos$BF.dB., xlab="SNP", ylab="BFaux (in dB)")
dev.off()