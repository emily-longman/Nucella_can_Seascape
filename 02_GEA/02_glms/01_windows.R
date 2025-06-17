# Generate windows for models

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
# install.packages(c('data.table', 'tidyverse', 'foreach', 'poolfstat', 'magrittr', 'reshape2', 'broom', 'SNPRelate'))
library(data.table)
library(tidyverse)
library(foreach)
library(poolfstat)
library(magrittr)
library(reshape2)
library(broom)
library(stats)

# Install and load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Load data

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# ================================================================================== #

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Define window and step size
win.bp = 5e4
step.bp = win.bp+1

# Generate windows (if a contig is less that the specified window size, the if statement, includes it)
wins <- foreach(chr.i=unique(snp.dt$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)

                  tmp <-  snp.dt %>%
                    filter(chr == chr.i)

                    if(max(tmp$pos) <= step.bp){
                      data.table(chr=chr.i, start=min(tmp$pos), end=max(tmp$pos))
                    } else if(max(tmp$pos) > step.bp){
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                } }

# Add column with window number
wins[,i:=1:dim(wins)[1]]

# ================================================================================== #

# Save windows
save(wins, file="data/processed/GEA/glms/glms_windows.RData")

