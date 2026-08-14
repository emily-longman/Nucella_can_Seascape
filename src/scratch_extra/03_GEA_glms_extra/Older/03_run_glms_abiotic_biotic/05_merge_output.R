# Analyze and graph output of model comparison

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
#install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)

# Prevent scientific notation
options(scipen=999)

# ================================================================================== #

# Load data

# Real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")

# ================================================================================== #

# Load Baypass input files for formatting

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
# Create SNP_id column
snp.meta <- snp.meta %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/output_chunks/', pattern = "GLM_output_chunk_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/output_chunks/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
o.all =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(o.all)

# Save merged data
save(o.all, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.output.all.Rdata"))

# ================================================================================== #

# Join real glm data with output
glm.o.tmp <- left_join(glm.real, o.all, by=c('chr', 'pos', 'SNP_id'))

# Order output
glm.o <- glm.o.tmp[order(snp.meta$SNP_id), ]
 
# ================================================================================== #

# Graph distribution of pval
pdf("output/figures/GEA/glms/Abiotic_biotic_pval_dist.pdf", width = 8, height = 8)
ggplot(glm.o, aes(x = p_val)) +
  geom_density() +
  theme_bw(base_size=25)
dev.off()

# Graph pval across genome
pdf("output/figures/GEA/glms/Abiotic_biotic_outliers.pdf", width = 14, height = 4)
ggplot(glm.o, aes(x = chr, y = -log(p_val))) +
  geom_point(alpha=0.6) +
  labs(x = "Position", y = "-log(P-val)") +
  theme_bw(base_size=25) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
dev.off()

