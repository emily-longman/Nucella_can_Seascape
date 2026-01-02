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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'doMC'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(doMC)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
win_group = as.numeric(args[1])

# ================================================================================== #

# State variable name
message(paste("Window group:", win_group))
# Load windows
load("data/processed/GEA/glms/glms_window_summary/windows.RData")
# Extract windows of interest
wins_group_i <- wins_guide_file_array %>% filter(.groups == win_group)

# Load ecological variables
Seascape_vars_names <- read.csv("guide_files/Seascape_vars_names.txt", header=F)
vars_subset <- Seascape_vars_names$V1[1:5]

# ================================================================================== #

# Register the multicore parallel backend
registerDoMC(20)

# Window summarization for each ecological variation
wins_sum <- foreach(var = vars_subset, .combine = "rbind", .errorhandling = "remove")%do%{

    # Load data
    load(paste0("data/processed/GEA/glms/glms_per_env_var/glm.collated_", var, ".Rdata"))

    # Create SNP_id column
    glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

    # Rank-normalize p-values
    glm.model.collated$rank <- rank(glm.model.collated$p_lrt)
    Lp <- length(glm.model.collated$p_lrt)
    glm.model.collated$rnp <- glm.model.collated$rank/Lp

    # Window summarization for each permutation and environmental var
    wins_sum <- foreach(perm.i=unique(glm.model.collated$perm),.combine="rbind", .errorhandling="remove")%dopar%{ 
    
        # State permutation number
        message(paste("Permutation #:", perm.i))

        # Filter data based on perm (0 = real data, 1 to 100 are permutations)
        tmp <- glm.model.collated %>% filter(perm == perm.i)

        # Start window summarization process
        win.out <- foreach(win.i=1:dim(wins_group_i)[1], .errorhandling = "remove", .combine = "rbind"
        )%dopar%{
    
        # State window
        message(paste("Window:", win.i, dim(wins_group_i)[1], sep = " / "))
  
        # Filter for a given window
        win_tmp <- tmp %>%
        filter(chr == wins_group_i[win.i]$chr) %>%
        filter(pos >= wins_group_i[win.i]$start & pos <= wins_group_i[win.i]$end)
    
        # P-values
        pr.i.0.01 <- c(0.01)
        pr.i.0.001 <- c(0.001)
    
        # Summarize for a given window
        win_tmp %>% 
            filter(!is.na(rnp)) %>%
            summarise(
              chr = wins_group_i[win.i]$chr,
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              perm = perm.i,
              variable = unique(variable),
              win = win.i,
              rnp.pr.0.01 = c(mean(rnp <= pr.i.0.01)),
              rnp.pr.0.001 = c(mean(rnp <= pr.i.0.001)),
              rnp.binom.p.0.01 = c(binom.test(sum(rnp <= pr.i.0.01), length(rnp), pr.i.0.01)$p.value),
              rnp.binom.p.0.001 = c(binom.test(sum(rnp <= pr.i.0.001), length(rnp), pr.i.0.001)$p.value),
              sum.rnp.0.01 = sum(rnp <= pr.i.0.01),
              sum.rnp.0.01 = sum(rnp <= pr.i.0.001),
              max.rnp = max(rnp),
              min.rnp = min(rnp),
              nSNPs = n()
            )
    }
}
}

# ================================================================================== #

# Generate folders and save output

# Folder name for chunk c
folder_name <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis")

# Save file for chunk w
file_name <- paste0("glm_window_chunks_", win_i)
save(wins_sum, file = paste0(folder_name, "/", file_name, ".Rdata") )

message("done")