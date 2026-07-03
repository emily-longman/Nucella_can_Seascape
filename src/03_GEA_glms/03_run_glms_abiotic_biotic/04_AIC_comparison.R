# Use AIC to calc p-val

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

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
w = as.numeric(args[1]) # Chunk

# Prevent scientific notation
options(scipen=999)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load data

# Real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")
# Perm data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.perm.Rdata")

# Load snp info with groups
load("data/processed/GEA/glms/snp.guide.file.Rdata")

# ================================================================================== #

# Subset real and perm data based on chunk

# Subset snp.guide. file
snp.chunk <- snp.guide.file %>% filter(.groups == w)
# Subset real data
glm.real.chunk <- glm.real %>% filter(SNP_id %in% snp.chunk$SNP_id)
# Subset perm data
glm.perm.chunk <- glm.perm %>% filter(SNP_id %in% snp.chunk$SNP_id)

# ================================================================================== #

# Identify which model is best for the real data

# Identify which column is the minimum
glm.real.chunk$minAIC <- names(glm.real.chunk[, 8:13])[apply(glm.real.chunk[, 8:13], MARGIN = 1, FUN = which.min)]
glm.real.chunk$minAIC_value <- do.call(pmin, c(glm.real.chunk[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.real.chunk <- glm.real.chunk %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# ================================================================================== #

# Identify which model is best for the perm data

# Identify which column is the minimum
glm.perm.chunk$minAIC_value <- do.call(pmin, c(glm.perm.chunk[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.perm.chunk <- glm.perm.chunk %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# ================================================================================== #

# Use permutations to calc p-val for SNPs where model with interaction is best fit 
o = foreach(i=1:dim(glm.real.chunk)[1], .combine = "rbind")%do%{
    
    # Extract real data for focal SNP
    real.tmp <- glm.real.chunk[i,]
    
    # Extract perm data for focal SNP
    perm.tmp <- glm.perm.chunk[which(glm.perm.chunk$SNP_id == glm.real.chunk$SNP_id[i]),]

    # Make data table and calculate p-val based on permutations (proportion of times that the perm data is less than the real data)
    data.frame(
          chr = unique(real.tmp$chr),
          pos = unique(real.tmp$pos),
          SNP_id = unique(real.tmp$SNP_id),
          real_minAIC = c(real.tmp$minAIC),
          real_AIC = c(real.tmp$deltaAIC),
          p_val = c(mean(perm.tmp$deltaAIC <= real.tmp$deltaAIC)))
}

# ================================================================================== #
# ================================================================================== #

# Use permutations to calc p-value

# Identify which SNPs in the real data have the model with the interaction as the best fit
#glm.real.chunk.int <- glm.real.chunk[which(glm.real.chunk$deltaAIC == 0),]
# Identify which SNPs in the real data dont have the model with the interaction as the best fit
#glm.real.chunk.no.int <- glm.real.chunk[which(glm.real.chunk$deltaAIC != 0),]

# Use permutations to calc p-val for SNPs where model with interaction is best fit 
#o.int = foreach(i=1:dim(glm.real.chunk.int)[1], .combine = "rbind")%do%{
#    
#    # Extract perm data for focal SNP
#    perm.tmp <- glm.perm.chunk[which(glm.perm.chunk$SNP_id == glm.real.chunk.int$SNP_id[i]),]
#
#    # Make data table and calculate p-val based on permutations
#    data.frame(
#          chr = unique(perm.tmp$chr),
#          pos = unique(perm.tmp$pos),
#          SNP_id = unique(perm.tmp$SNP_id),
#          p_val = c(1-length(which(perm.tmp$deltaAIC == 0))/dim(perm.tmp)[1]))
#}

# For all models where the interaction is not the best fit, set p-val equal to 1
#o.no.int <- glm.real.chunk.no.int[,1:3] %>% mutate(p_val = 1)

# Join data
#o <- rbind(o.int, o.no.int)

# ================================================================================== #

# Generate folders and save output

# Folder name
folder_name <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/output_chunks")
if (!dir.exists(folder_name)) {dir.create(folder_name)}

# Save file for chunk w
file_name <- paste("GLM_output_chunk_", w, sep = "")
save(o, file = paste(folder_name, "/" , file_name, ".Rdata", sep = "") )

message("done")

# ================================================================================== #
# ================================================================================== #

# Test
#realtmp <- glm.model.collated.real[6,]
#permtmp <- glm.model.collated.perm1.50[which(glm.model.collated.perm1.50$SNP_id==realtmp$SNP_id),]

#1-length(which(permtmp$minAIC == "AIC_dem_both_int"))/50

#pdf("output/figures/GEA/glms/test.pdf", width = 14, height = 6)
#ggplot(permtmp, aes(x = deltaAIC)) +
#    geom_density(alpha = 0.7, lwd = 1) +
#    theme_bw(base_size=30)
#dev.off()