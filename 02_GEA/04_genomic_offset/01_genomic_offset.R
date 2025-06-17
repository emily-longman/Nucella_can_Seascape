# Perform genomic offset analyses with Bio-Oracle data
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'terra', 'raster', 'SeqArray'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(terra)
library(raster)

# Load gradientForest
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
library(gradientForest)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Load poolobject
load("data/processed/pop_structure/pooldata.RData")

# Extract data for SNPs of interest
str(pooldata)
head(pooldata)


# ================================================================================== #


# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Read in significant SNPs
baypass_pos_bonf_sig_SNPs <- read.table("data/processed/GEA/baypass/abiotic/baypass_pos_bonf_sig_SNPs.txt", header=T)

# Extract SNP metadata from genofile
genofile_SNP <- data.table(
  variant.id = seqGetData(genofile, "variant.id"),
  chr = seqGetData(genofile, "chromosome"),
  pos = seqGetData(genofile, "position"))

# Filter genofile to only significant SNPs
sig_SNPs <- genofile_SNP %>%
  semi_join(baypass_pos_bonf_sig_SNPs, by = c("chr", "pos"))

# Set filter in GDS to only those SNPs
seqSetFilter(genofile, variant.id = sig_SNPs$variant.id)







#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")


# Load SNPs of interest (baypass rnp outlier SNPs: p=0.001)
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.001.RData")
# Create SNP_id column for outlier SNP list
SNPs.Interest.pval.001 <- SNPs.Interest.pval.001 %>% mutate(SNP_id = paste(chr, pos, sep = "_"))


# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        id=seqGetData(genofile, "variant.id")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# Filter for only outlier SNPs (NOTE: somehow getting only 61,896 SNPs)
sig_snps <- snp.dt %>% filter(snp.dt$SNP_id %in% SNPs.Interest.pval.001$SNP_id)

#--------------------------------------------------------------------------------

# Set filter in GDS to only those SNPs
seqSetFilter(genofile, variant.id = sig_snps$id)

seqGetData(genofile, "genotype") -> genotypes_i

seqGetData(genofile, "sample.id") -> sampsids_raw

