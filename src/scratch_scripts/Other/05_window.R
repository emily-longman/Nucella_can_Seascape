# Generate windows for baypass outlier analyses.

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
#install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
#install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen'))
library(data.table)
library(dplyr)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Get files from bash script
argv <- commandArgs(T)
XtX_path<-argv[1]
snp.meta_path<-argv[2]

# Read in Baypass results
XtX <- read.table(XtX_path, header=T)
snp.meta <- read.table(snp.meta_path, header=F)


# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Generate windows
wins <- foreach(chr.i=unique(SNP.XtX.dt$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  tmp <-  inner.rnf %>%
                    filter(chr == chr.i)

                    if(max(tmp$pos) <= step.bp){
                      data.table(chr=chr.i,
                             start=min(tmp$pos),
                             end=max(tmp$pos))  
                    } else if(max(tmp$pos) > step.bp){
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                } }

wins[,i:=1:dim(wins)[1]]

# Then do an alpha of 0.01 or 0.05

save.image("data/processed/GEA/baypass/Baypass_windows.RData")