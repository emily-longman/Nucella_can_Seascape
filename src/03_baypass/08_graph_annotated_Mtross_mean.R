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
out_fig_dir <- paste("output/figures/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

# Load baypass BF data (bf.Mtross.mean.sum)
load("data/processed/baypass/biotic/bf.Mtross.mean.sum.Rdata")

# Load annotated outliers
bf.Mtross.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.Mtross.mean.sum.outliers.annotated.csv", header=T)

# ================================================================================== #

# How many unique genes - 2237
length(unique(bf.Mtross.mean.sum.outliers.annotated$Gene_Name))

# Join annotated outliers to full dataset
Mtross <- left_join(bf.Mtross.mean.sum, bf.Mtross.mean.sum.outliers.annotated, 
    by = c("chr", "pos", "allele1", "allele2", "MRK", "bf_db.mean", "bf_db.median", "bf_db.var", "eBPis.mean", "eBPis.median", "eBPis.var"))

# Color points in gene g26813
Mtross$color <- ifelse(Mtross$Gene_Name == "g26813", "g26813", "not.g26813")


# Graph BF with 0.001 POD threshold - only BF > 0
pdf("output/figures/baypass/baypass_BF_Mtross_mean_repmeans_posBF.pdf", width = 12, height = 8)
ggplot(bf.Mtross.mean.sum[which(bf.Mtross.mean.sum$bf_db.mean>0),], aes(y=bf_db.mean, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))
dev.off()