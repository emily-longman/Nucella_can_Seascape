# Summarize glms output and compare across the environmental variables

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
#install.packages(c('data.table', 'tidyverse', 'ggplot2', 'foreach', 'dplyr'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Get names of enviro variables
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Extract just names
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# ================================================================================== #

# Model enrichment

# Summarize the permutation data across the environmental variables
all_ratios <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
    # State variable name
    message(i)

    # Load data
    load(paste0("data/processed/GEA/glms/glms_per_env_var/glm.collated_", i, ".Rdata") )

    # Graph pval distribution - abiotic
    #pdf(paste0("output/figures/GEA/glms/glm_pval_dist_abiotic_", i, ".pdf"), width = 8, height = 8)
    #ggplot(glm.model.collated, aes(x=p_lrt, group=factor(perm), color=factor(perm))) + geom_density() +
    #scale_color_manual(values = c("red", rep("grey", 100))) +
    #xlab("GLM P-values") + ylab("Number of SNPs") +
    #theme_bw() + theme(legend.position = "none")
    #dev.off()

    # Number of permutations
    n_perm = length(unique(glm.model.collated$perm[which(glm.model.collated$perm>0)]))
    message(n_perm)

    # Check if the environmental model beats the demographic model by comparing AIC
    glm.model.collated <- glm.model.collated %>% mutate(AICdiff = if_else(AIC_dem_env < AIC_dem, 1, 0)) 

    # Calculate relative rate (rr) of model enrichment (i.e., the number of SNPs where the env dem model was found as the better model)
    aic_sum <- glm.model.collated %>% group_by(perm) %>% summarise(rr = sum(AICdiff == 1, na.rm = TRUE), .groups = "drop")

    # Separate by real and perm
    real_data <- aic_sum %>% filter(perm == 0) %>% rename(rr_real = rr)
    perm_data <- aic_sum %>% filter(perm > 0) %>% rename(rr_perm = rr)

    # rr of real data
    rr_real <- real_data$rr_real

    # Join the datasets and compare the rr between real and perm
    ratios <- perm_data %>%
    mutate(rr_ratio = (rr_real / rr_perm), rr_ratio_log2 = log2(rr_real / rr_perm), variable = i)
    
    # Return ratios
    return(ratios)

    # Remove glm.model.collated for that env var
    rm(glm.model.collated)
}

# Write table
write.csv(all_ratios, "data/processed/GEA/glms/glms_summary/Bio-oracle.rr.csv", row.names=FALSE)

