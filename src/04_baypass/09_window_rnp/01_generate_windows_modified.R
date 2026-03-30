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
out_dir <- paste("data/processed/baypass/window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
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

# How many SNPs are on each contig:
SNPS_density <- snp.dt %>% group_by(chr) %>% summarize(n=n())

# Graph SNP density
pdf("output/figures/baypass/window_summary/glm_chr_nSNP_density.pdf", width = 8, height = 8)
ggplot(SNPS_density, aes(x=n))+ geom_density() + xlim(0,150)
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Extract SNP info for all SNPs
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  pooldata.snp.info

# Rename columns
names(pooldata.snp.info)[1:2] = c("chr","pos")

# Make snp_id column
pooldata.snp.info<- pooldata.snp.info %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

# Filter GDS snp.dt to only sites in pooldata snp.info
snp.dt.filt <- snp.dt %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)

######

# How many SNPs are on each contig:
SNPS_density_filt <- snp.dt.filt %>% group_by(chr) %>% summarize(n=n())

# Graph SNP density
pdf("output/figures/baypass/window_summary/glm_chr_nSNP_filt_density.pdf", width = 8, height = 8)
ggplot(SNPS_density_filt, aes(x=n))+ geom_density() + xlim(0,150)
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# ================================================================================== #

# Read in SNP data
snpdet <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snpdet) <- c("chr", "pos", "allele1", "allele2")


SNPS_snpdet_sum <- snpdet %>% group_by(chr) %>% summarize(n=n())

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 100000
step.bp <- 50000

# Generate windows
wins <- foreach(chr.i=unique(snp.dt.filt$chr), .combine="rbind", .errorhandling="remove")%do%{
      # State chromosome
      message(chr.i)

      # Filter data for focal chromosome
      tmp <- snp.dt.filt %>%
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
write.table(wins_guide_file_array, "guide_files/wins_guide_file_array.txt", col.names = F, row.names = F, quote = F)
# Note guide_file_array has dimensions: 12318, 5, - 493 groups

# ================================================================================== #

# Save windows
save(wins_guide_file_array, file="data/processed/baypass/window_summary/windows.RData")