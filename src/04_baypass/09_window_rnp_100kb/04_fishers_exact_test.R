# Analyze windows

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'doMC'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# Load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
outlier.win.SNPs.annotated.ph_mean <- read.csv("data/processed/baypass/window_summary/outlier.win.SNPs.annotated.ph_mean.csv", header=T)
outlier.win.SNPs.annotated.Mtross_mean <- read.csv("data/processed/baypass/window_summary/outlier.win.SNPs.annotated.Mtross_mean.csv", header=T)
outlier.win.SNPs.annotated.McaliIntThk <- read.csv("data/processed/baypass/window_summary/outlier.win.SNPs.annotated.McaliIntThk.csv", header=T)

# Summarize each list based on annotation
ph_ann_sum <- outlier.win.SNPs.annotated.ph_mean %>% count(Annotation) %>% rename(ph = n)
Mtross_ann_sum <- outlier.win.SNPs.annotated.Mtross_mean %>% count(Annotation) %>% rename(Mtross = n)
McaliIntThk_ann_sum <- outlier.win.SNPs.annotated.McaliIntThk %>% count(Annotation) %>% rename(McaliIntThk = n)

# Join data
ann <- full_join(ph_ann_sum, Mtross_ann_sum)
ann <- full_join(ann, McaliIntThk_ann_sum)
ann[is.na(ann)] <- 0
row.names(ann) <- ann$Annotation
ann <- ann[,c(2:4)]

pdf("output/figures/baypass/window_summary/mosaic_plot.pdf", width = 12, height = 6)
mosaicplot(ann, color = TRUE)
dev.off()

# Fishers exact test (Monte Carlo simulation with 2,000 simulations)
ftest <- fisher.test(ann, simulate.p.value = TRUE, B = 2000)

# Chi-sq test
chisq.test(ann)

# So they are different, but there are a ton of categories so that feels less surprising