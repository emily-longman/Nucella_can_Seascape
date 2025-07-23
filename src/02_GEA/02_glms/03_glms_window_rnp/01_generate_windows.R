# Create windows for GLM summarization

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr'))
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

# Load GLM data
load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 1e5
step.bp <- 5e8

# How many SNPs are on each contig:
SNPS_density <- glm.model.collated %>% group_by(chr) %>% summarize(n=n())
# Graph
pdf("output/figures/GEA/glms/glms_window_summary/glm_pval_density.pdf", width = 8, height = 8)
ggplot(SNPS_density, aes(x=n))+ geom_density()
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# Generate windows (note: only windows with the number of SNPs in that window >= 5)
wins <- foreach(chr.i=unique(glm.model.collated$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  
                  tmp <- glm.model.collated %>%
                    filter(chr == chr.i)
                  
                  nSNPs=dim(tmp)[1]
                  
                  if(nSNPs >= 5){
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

# Check dimensions
dim(wins)

# ================================================================================== #

# Save windows
save(wins, file="data/processed/GEA/glms/glms_window_summary/windows.RData")
