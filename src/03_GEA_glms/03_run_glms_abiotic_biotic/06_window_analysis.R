# Window analysis of abiotic biotic glms

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
load("data/processed/baypass/window_summary/windows_100kb.RData")
# Extract windows of interest based on win_group
wins_group_i <- wins_guide_file_array %>% filter(.groups == win_group)

# Load real data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.real.Rdata")

# Load output
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.output.all.Rdata")

# ================================================================================== #

# Load Baypass input files for formatting

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
# Create SNP_id column
snp.meta <- snp.meta %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Join real glm data with output
glm.o.tmp <- left_join(glm.real, o.all, by=c('chr', 'pos', 'SNP_id'))

# Order output
glm.o <- glm.o.tmp[order(snp.meta$SNP_id), ]

# ================================================================================== #

# Rank-normalize bf      (note: high bf should be associated with a low rank)
glm.o$rank <- rank(glm.o$p_val)
Lp <- length(glm.o$p_val)
glm.o$rnp <- glm.o$rank/Lp

# ================================================================================== #

# Window summarization

# Start window summarization process for windows within window group
win.out <- foreach(window.w=1:dim(wins_group_i)[1], .combine = "rbind", .errorhandling = "remove")%do%{

        # State window
        message(paste0("Window ", window.w, " / ", dim(wins_group_i)[1], " for window group"))
        message(paste0("i.e., Window #: ", wins_group_i$i[window.w]))
  
        # Filter glm data for a given window
        win_tmp <- glm.o %>%
        filter(chr == wins_group_i$chr[window.w]) %>%
        filter(pos >= wins_group_i$start[window.w] & pos <= wins_group_i$end[window.w])

        # Set p
        pr.i = 0.05

        # Summarize for a given window
        win_tmp %>% 
            filter(!is.na(rnp)) %>%
            summarise(
              chr = unique(chr),
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              p_val_mean = mean(p_val),
              p_val_min = min(p_val),
              p_val_max = max(p_val),
              win = wins_group_i$i[window.w],
              rnp.pr = c(mean(rnp <= pr.i)),
              rnp.binom.p = c(binom.test(sum(rnp <= pr.i), length(rnp), pr.i)$p.value),
              sum.rnp = sum(rnp <= pr.i),
              max.rnp = max(rnp),
              min.rnp = min(rnp),
              nSNPs = n()
            )
}

# ================================================================================== #

# Generate folders and save output

# Folder name for window group i
folder_name <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic")

# Save file for chunk w
file_name <- paste0("glm_window_chunks_", win_group)
save(wins_sum, file = paste0(folder_name, "/", file_name, ".Rdata") )

message("done")