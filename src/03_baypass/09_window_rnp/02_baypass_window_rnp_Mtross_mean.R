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

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
win_group = as.numeric(args[1])

# ================================================================================== #

# State variable name
message(paste("Window group:", win_group))
# Load windows
load("data/processed/baypass/window_summary/windows.RData")
# Extract windows of interest based on win_group
wins_group_i <- wins_guide_file_array %>% filter(.groups == win_group)

# Load mean bf data from 5 baypass runs
load("data/processed/baypass/biotic/bf.Mtross.mean.sum.Rdata")

# Load POD thresholds
load("data/processed/baypass/biotic/Mtross_mean_POD_thr.Rdata")

# ================================================================================== #

# Rank-normalize bf (note: high bf should be associated with a low rank)
bf.Mtross.mean.sum$rank <- rank(-bf.Mtross.mean.sum$bf_db.mean)
Lp <- length(bf.Mtross.mean.sum$bf_db.mean)
bf.Mtross.mean.sum$rnp <- bf.Mtross.mean.sum$rank/Lp

# ================================================================================== #

# Use the POD threshold to come up with p-val

# Create the ECDF (empirical cumulative distribution functio) function
my_ecdf <- ecdf(bf.Mtross.mean.sum$bf_db.mean)

# Find the probability for a given value
probability <- my_ecdf(bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)])

# ================================================================================== #

# Window summarization

# Start window summarization process for windows within window group
win.out <- foreach(window.w=1:dim(wins_group_i)[1], .combine = "rbind", .errorhandling = "remove")%do%{

        # State window
        message(paste0("Window ", window.w, " / ", dim(wins_group_i)[1], " for window group"))
        message(paste0("i.e., Window #: ", wins_group_i$i[window.w]))
  
        # Filter glm data for a given window
        win_tmp <- bf.Mtross.mean.sum %>%
        filter(chr == wins_group_i$chr[window.w]) %>%
        filter(pos >= wins_group_i$start[window.w] & pos <= wins_group_i$end[window.w])
    
        # Summarize for a given window
        win_tmp %>% 
            filter(!is.na(rnp)) %>%
            summarise(
              chr = unique(chr),
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              bf_mean = mean(bf_db.mean),
              bf_min = min(bf_db.mean),
              bf_max = max(bf_db.mean),
              win = wins_group_i$i[window.w],
              rnp.POD = c(mean(rnp <= probability)),
              rnp.binom.POD = c(binom.test(sum(rnp <= probability), length(rnp), probability)$p.value),
              sum.rnp.POD = sum(rnp <= probability),
              max.rnp = max(rnp),
              min.rnp = min(rnp),
              nSNPs = n()
            )
}

# ================================================================================== #

# Generate folders and save output

# Folder name for window group i
folder_name <- paste("data/processed/baypass/window_summary/window_chunk_analysis_Mtross_mean")

# Save file for chunk w
file_name <- paste0("window_chunks_", win_group)
save(win.out, file = paste0(folder_name, "/", file_name, ".Rdata") )

message("done")