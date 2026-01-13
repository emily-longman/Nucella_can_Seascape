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
#install.packages(c('data.table', 'tidyverse', 'ggplot2', 'foreach', 'dplyr'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(foreach)
library(dplyr)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
env_var = as.character(args[1]) #Environmental variable

# ================================================================================== #

# Generate output directories
out_dir <- paste("/gpfs2/scratch/elongman/Nucella_can_Seascape/data/processed/GEA/glms/glms_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Model enrichment - Summarize the permutation data across the environmental variables

# State variable name
message(env_var)

# Load data
load(paste0("/gpfs3/scratch/elongman/glms_per_env_var/glm.collated_", env_var, ".Rdata") )

# Bin data
hist.obj.env = foreach(i = 0:max(glm.model.collated$perm), .combine = "rbind")%do%{
hist(glm.model.collated$p_lrt[glm.model.collated$perm == i], breaks = 100) -> hist.obj
data.frame(
  hist.obj$mids,
  hist.obj$counts,
  perm = i) -> o
}

# Graph pval distribution
pdf(paste0("output/figures/GEA/glms/glm_pval_dist_log_scale_", env_var, ".pdf"), width = 8, height = 8)
ggplot(hist.obj.env, aes(x=(hist.obj.mids),y=hist.obj.counts, group=perm, color=perm==0)) +
  geom_line(aes(alpha=perm==0 ), linewidth = 2.5) + 
  scale_alpha_manual(values = c(0.1, 1)) +
  scale_size_manual(values = c(0.7, 1.3)) +
  scale_color_manual(values = c("red","black")) +
  labs(title = paste0(env_var, " P-value distribution"), x = "GLM P-values", y = "Number of SNPs") +
  theme_bw(base_size = 24) +  theme(legend.position = "none") +
  scale_x_log10()
dev.off()

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
mutate(rr_ratio = (rr_real / rr_perm), rr_ratio_log2 = log2(rr_real / rr_perm), variable = env_var)

# ================================================================================== #

# Write table
write.csv(ratios, file = paste("/gpfs2/scratch/elongman/Nucella_can_Seascape/data/processed/GEA/glms/glms_summary/Vars_rr_", env_var, ".csv", sep = ""), row.names=FALSE)
