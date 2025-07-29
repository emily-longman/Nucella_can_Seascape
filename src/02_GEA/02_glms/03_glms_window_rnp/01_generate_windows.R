# Create windows for GLM summarization

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
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/GEA/glms/glms_window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load SNP data (3,095 outlier SNPs)
baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.01", header=T) 

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 1e5 #100,000
step.bp <- 5e4 #50,000

# How many SNPs are on each contig:
SNPS_density <- baypass_POD_sig_SNPs %>% group_by(chr) %>% summarize(n=n())

# Graph
pdf("output/figures/GEA/glms/glms_window_summary/glm_pval_density.pdf", width = 8, height = 8)
ggplot(SNPS_density, aes(x=n))+ geom_density()
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# Generate windows (note: only windows with the number of SNPs in that window >= 2)
wins <- foreach(chr.i=unique(baypass_POD_sig_SNPs$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  
                  tmp <- baypass_POD_sig_SNPs %>%
                    filter(chr == chr.i)
                  
                  nSNPs=dim(tmp)[1]
                  
                  if(nSNPs >= 2){
                    o =
                      data.table(chr=chr.i,
                                 nSNPs=dim(tmp)[1],
                                 start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                                 end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                    return(o)
                    
                  }   
                  else {message("fails nSNPs filter")}
                }

# Add window index
wins[,i:=1:dim(wins)[1]]

# Check dimensions - 335 windows
dim(wins)

# ================================================================================== #

# Save windows
save(wins, file="data/processed/GEA/glms/glms_window_summary/windows_SNPs_threshold_0.01.RData")
