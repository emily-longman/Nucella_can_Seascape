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

# Register the multicore parallel backend
registerDoMC(10)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
win_group = as.numeric(args[1])

# ================================================================================== #

# State variable name
message(paste("Window group:", win_group))
# Load windows
load("data/processed/GEA/glms/glms_window_summary/windows.RData")
# Extract 13 windows of interest based on win_group
wins_group_i <- wins_guide_file_array %>% filter(.groups == win_group)

# Load ecological variables
#Seascape_vars_names <- read.csv("guide_files/Seascape_vars_names.txt", header=F)
#vars_subset <- Seascape_vars_names$V1[1:5]

# Load data
load("data/processed/GEA/glms/glms_per_env_var/glm.collated_ph_mean.Rdata")

# ================================================================================== #

# Window summarization for ph mean

# Window summarization for each permutation and environmental var
#wins_sum <- foreach(perm.i=unique(glm.model.collated$perm),.combine="rbind", .errorhandling="remove")%dopar%{ 
    perm.i=0
    # State permutation number
    message(paste("Permutation #:", perm.i))

    # Filter glm data based on perm (0 = real data, 1 to 100 are permutations)
    tmp <- glm.model.collated %>% filter(perm == perm.i)

        # Rank-normalize p-values
        tmp$rank <- rank(tmp$p_lrt)
        Lp <- length(tmp$p_lrt)
        tmp$rnp <- tmp$rank/Lp

        # Start window summarization process for windows within window group
        win.out <- foreach(window.w=1:dim(wins_group_i)[1], .combine = "rbind", .errorhandling = "remove")%dopar%{
    
        # State window
        message(paste0("Window: ", window.w, " / ", dim(wins_group_i)[1]))
        message(paste0("i.e., Window #: ", wins_group_i$i[window.w]))
  
        # Filter glm data for a given window
        win_tmp <- tmp %>%
        filter(chr == wins_group_i$chr[window.w]) %>%
        filter(pos >= wins_group_i$start[window.w] & pos <= wins_group_i$end[window.w])
    
        # P-values
        pr.i.0.01 <- c(0.01)
        pr.i.0.001 <- c(0.001)
    
        # Summarize for a given window
        win_tmp %>% 
            filter(!is.na(rnp)) %>%
            summarise(
              chr = unique(chr),
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              variable = unique(variable),
              perm = perm.i,
              win = wins_group_i$i[window.w],
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
    #return(win.out)
#}

# ================================================================================== #

# Generate folders and save output

# Folder name for window group i
folder_name <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_ph_mean_real")

# Save file for window group
file_name <- paste0("glm_window_chunks_", win_group)
save(win.out, file = paste0(folder_name, "/", file_name, ".Rdata") )

message("done")