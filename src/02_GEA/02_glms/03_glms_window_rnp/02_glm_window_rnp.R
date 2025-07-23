# Create windows 

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/GEA/glms/glms_window_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load GLM data
load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")
#load("data/processed/GEA/glms/glms_output/glm.model.collated.test.Rdata")

# Load windows
load("data/processed/GEA/glms/glms_window_summary/windows.RData")

# ================================================================================== #

# Format data

# Create SNP_id column
glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Rank-normalize p-values

# Filter for just the real data
real_data <- glm.model.collated %>% filter(perm == 0)

# Rank pvalues
real_data_rank <- foreach(var.i=unique(real_data$variable),.combine="rbind", .errorhandling="remove")%do%{ 
        # Filter for just one environmental var
        tmp <- real_data %>% filter(variable == var.i)
        # Rank p values
        tmp %>% mutate(
            rank=rank(p_lrt),
            Lp = length(p_lrt),
            rnp=rank/Lp)
                }


# ================================================================================== #

# Window summarization for each environmental var
real_data_win_sum <- foreach(var.i=unique(real_data$variable),.combine="rbind", .errorhandling="remove")%do%{ 
    
    # Filter for just one environmental var
    tmp <- real_data %>% filter(variable == var.i)

    # Start window summarization process
    win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
    )%do%{
  
    message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
    win.tmp <- real_data_rank %>%
        filter(chr == wins[win.i]$chr) %>%
        filter(pos >= wins[win.i]$start & pos <= wins[win.i]$end)
  
    pr.i.0.01 = 0.01
    pr.i.0.001 = 0.001
  
    win.tmp %>% 
        filter(!is.na(rnp)) %>%
        summarise(
              chr = wins[win.i]$chr,
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              variable = var.i,
              win = win.i,
              rnp.pr.0.01 = c(mean(rnp <= pr.i.0.01)),
              rnp.pr.0.001 = c(mean(rnp <= pr.i.0.001)),
              rnp.binom.p.0.01 = c(binom.test(sum(rnp <= pr.i.0.01), length(rnp), pr.i.0.01)$p.value),
              rnp.binom.p.0.001 = c(binom.test(sum(rnp <= pr.i.0.001), length(rnp), pr.i.0.001)$p.value),
              sum.rnp.0.01 = sum(rnp <= pr.i.0.01),
              sum.rnp.0.001 = sum(rnp <= pr.i.0.001),
              max.p_lrt = max(p_lrt),
              min.rnp = min(rnp),
              nSNPs = n()
            )  -> win.out
    }
}

# ================================================================================== #

# Save glm rnp for windows
save(real_data_win_sum, file="data/processed/GEA/glms/glms_window_summary/glm_windows_output.RData")
