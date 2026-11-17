# Create windows

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)

# Install and load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/baypass/window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Read in SNP data
snpdet <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snpdet) <- c("chr", "pos", "allele1", "allele2")

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 100000
step.bp <- 50000

# Generate windows
wins <- foreach(chr.i=unique(snpdet$chr), .combine="rbind", .errorhandling="remove")%do%{
      # State chromosome
      message(chr.i)

      # Filter data for focal chromosome
      tmp <- snpdet %>%
      filter(chr == chr.i)

      # Number of SNPs on chromosome
      nSNPs=dim(tmp)[1]

      # For only chromosomes with >= 100
      if(nSNPs >= 100){
      o = data.table(
        chr=chr.i,
        nSNPs=dim(tmp)[1],
        start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
        end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)

      # Return output
      return(o)

      }   
      else {message("fails nSNPs filter")}
}

# Add window index
wins[,i:=1:dim(wins)[1]]

# Check dimensions - 12,318 windows
dim(wins)

# ================================================================================== #

# Group the backbone names into ~500
group(wins, n=25, method = "greedy") -> wins_guide_file_array

# Write the table
write.table(wins_guide_file_array, "guide_files/wins_10kb_guide_file_array.txt", col.names = F, row.names = F, quote = F)
# Note guide_file_array has dimensions: 12318, 5, - 493 groups

# ================================================================================== #

# Save windows
save(wins_guide_file_array, file="data/processed/baypass/window_summary/windows_100kb.RData")