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

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load data

# Real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.real.Rdata")
glm.real <- glm.model.collated.real
# Perm data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.perm.Rdata")
glm.perm <- perm

# ================================================================================== #

# Identify which model is best for the real data

# Identify which column is the minimum
#glm.real$minAIC <- names(glm.real[, 8:13])[apply(glm.real[, 8:13], MARGIN = 1, FUN = which.min)]
glm.real$minAIC_value <- do.call(pmin, c(glm.real[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.real <- glm.real %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# ================================================================================== #

# Identify which model is best for the perm data

# Identify which column is the minimum
glm.perm$minAIC_value <- do.call(pmin, c(glm.perm[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.perm <- glm.perm %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# ================================================================================== #

# Use permutations to calc p-value

# Identify which SNPs in the real data have the model with the interaction as the best fit
glm.real.int <- glm.real[which(glm.real$deltaAIC == 0),]
# Identify which SNPs in the real data dont have the model with the interaction as the best fit
glm.real.no.int <- glm.real[which(glm.real$deltaAIC != 0),]

# Use permutations to calc p-val for SNPs where model with interaction is best fit 
o.int = foreach(i=1:dim(glm.real.int)[1], .combine = "rbind")%do%{
    
    # Extract perm data for focal SNP
    perm.tmp <- glm.perm[which(glm.perm$SNP_id == glm.real.int$SNP_id[i]),]

    # Make data table and calculate p-val based on permutations
    data.frame(
          chr = unique(perm.tmp$chr),
          pos = unique(perm.tmp$pos),
          SNP_id = unique(perm.tmp$SNP_id),
          p_val = c(1-length(which(perm.tmp$deltaAIC == 0))/dim(perm.tmp)[1]))
}

# For all models where the interaction is not the best fit, set p-val equal to 1
o.no.int <- no.int[,1:3] %>% mutate(p_val = 1)

# Join data
o <- rbind(o.int, o.no.int)

# ================================================================================== #

# Save output
save(o, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.output.Rdata"))

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