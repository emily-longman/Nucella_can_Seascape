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
out_dir <- paste("data/processed/GEA/glms/glms_window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Extract SNP data from GDS (14,897,468 SNPs)
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# Load SNP data (3,095 outlier SNPs)
#baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.01", header=T) 

# ================================================================================== #

# How many SNPs are on each contig:
SNPS_density <- snp.dt %>% group_by(chr) %>% summarize(n=n())

# Graph SNP density
pdf("output/figures/GEA/glms/glms_window_summary/glm_pval_density.pdf", width = 8, height = 8)
ggplot(SNPS_density, aes(x=n))+ geom_density() + xlim(0,150)
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 100000
step.bp <- 50000

# Generate windows (note: only chromosomes with the number of SNPs in that window >= 5)
wins <- foreach(chr.i=unique(snp.dt$chr), .combine="rbind", .errorhandling="remove")%do%{
      # State chromosome
      message(chr.i)

      # Filter data for focal chromosome
      tmp <- snp.dt %>%
      filter(chr == chr.i)

      # Number of SNPs on chromosome
      nSNPs=dim(tmp)[1]

      # For only chromosomes with >= 5
      if(nSNPs >= 5){
      o = data.table(chr=chr.i,
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

# Check dimensions - 12928 windows
dim(wins)

# ================================================================================== #

# Group the backbone names into 1000 (Note: 12928 wins /1000 array = 12.9)
group(wins, n=13, method = "greedy") -> wins_guide_file_array

# Write the table
write.table(wins_guide_file_array, "guide_files/wins_guide_file_array.txt", col.names = F, row.names = F, quote = F)
# Note guide_file_array has dimensions:  12928, 6 - 995 groups

# ================================================================================== #

# Save windows
save(wins_guide_file_array, file="data/processed/GEA/glms/glms_window_summary/windows.RData")