# Summarize glms output and compare across the environmental variables

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
#install.packages(c('rprojroot'))
#library(rprojroot)
# Specify root path
#root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
#setwd(root_path)

# ================================================================================== #

# Load packages
#install.packages(c('data.table', 'tidyverse', 'ggplot2', 'foreach', 'doMC', 'dplyr'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(foreach)
library(doMC)
library(dplyr)

# ================================================================================== #

# Generate output directories
out_dir <- paste("/gpfs2/scratch/elongman/Nucella_can_Seascape/data/processed/GEA/glms/glms_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Get names of enviro variables
vars_names <- read.csv("guide_files/Seascape_vars_names.txt", header=F)

# Extract just names
vars_names <- vars_names$V1

# ================================================================================== #

# Model enrichment
registerDoMC(20)

# Summarize the permutation data across the environmental variables
all_ratios <- foreach(i=vars_names, .combine="rbind")%dopar%{
    # State variable name
    message(i)

    # Load data
    load(paste0("/gpfs3/scratch/elongman/glms_per_env_var/glm.collated_", i, ".Rdata") )

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


# ================================================================================== #

# Write table
write.csv(all_ratios, "/gpfs2/scratch/elongman/Nucella_can_Seascape/data/processed/GEA/glms/glms_summary/Abiotic_biotic_vars_rr.csv", row.names=FALSE)

