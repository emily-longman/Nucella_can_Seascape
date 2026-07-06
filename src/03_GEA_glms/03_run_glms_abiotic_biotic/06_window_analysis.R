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

# Rank-normalize
#glm.o$rank <- rank(glm.o$p_val)
#Lp <- length(glm.o$p_val)
#glm.o$rnp <- glm.o$rank/Lp
# Load and average xtx

# ================================================================================== #

# Did once when testing and then loaded for subsequent array

# NOTE these xtx values are for all 19 sites

# Read in SNP data
#snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
#colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load xtx output for 5 replicate Baypass runs
#baypass.xtx <- foreach(i=1:5, .combine = rbind)%do%{
#    message(i)
#    tmp <- fread(paste("data/processed/baypass/xtx/NC_run", i, "_summary_pi_xtx.out", sep=""))
#    tmp[,rep:=i]
#    tmp <- cbind(snp.meta, tmp)
#    return(tmp)
#}

# Rename p val
#baypass.xtx <- baypass.xtx %>% rename(log10.1.pval. = "log10(1/pval)")

# Average across replicate runs
#baypass.xtx <- baypass.xtx %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
#    reframe(M_P_mean = mean(M_P), SD_P_mean = mean(SD_P), M_XtX_mean = mean(M_XtX), 
#    SD_XtX_mean = mean(SD_XtX), XtXst_mean = mean(XtXst), log10.1.pval_mean = mean(log10.1.pval.))

# Save
#save(baypass.xtx, file = "data/processed/baypass/baypass.xtx.sum.RData")
load("data/processed/baypass/baypass.xtx.sum.RData")

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

        # Extract xtx data for a given window
        win_xtx <- baypass.xtx %>%
        filter(chr == wins_group_i$chr[window.w]) %>%
        filter(pos >= wins_group_i$start[window.w] & pos <= wins_group_i$end[window.w])

        # Set p
        pr.i = 0.05
        #pr.i = 0.01
  
  
        # Create matrix for significance of focal annotation
        win_tmp_sum <- matrix(c(
          length(which(win_tmp$p_val < 0.05 & win_tmp$real_deltaAIC == 0)),
          length(which(win_tmp$p_val >= 0.05 & win_tmp$real_deltaAIC == 0)),
          length(which(win_tmp$p_val < 0.05 & win_tmp$real_deltaAIC != 0)),
          length(which(win_tmp$p_val >= 0.05 & win_tmp$real_deltaAIC != 0))),
          nrow = 2)

        # Fishers exact test
        ftest_w <- fisher.test(win_tmp_sum)

        # Create data frame with output (i, p.value, odds ratio, and 95% CI)
        data.frame(
          chr = unique(win_tmp$chr),
          pos_mean = mean(win_tmp$pos),
          pos_min = min(win_tmp$pos),
          pos_max = max(win_tmp$pos),
          win = wins_group_i$i[window.w],
          n_intAIC_bestsig = length(which(win_tmp$p_val < 0.05 & win_tmp$real_deltaAIC == 0)),
          p.fet = ftest_w$p.value,
          OR = ftest_w$estimate,
          lci = ftest_w$conf.int[1],
          uci = ftest_w$conf.int[2],
          meanxtx = mean(win_xtx$XtXst_mean),
          nSNPs = dim(win_tmp)[1]
        )

}

# ================================================================================== #

# Generate folders and save output

# Folder name for window group i
folder_name <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic_FET")

# Save file for chunk w
file_name <- paste0("glm_window_chunks_", win_group)
save(win.out, file = paste0(folder_name, "/", file_name, ".Rdata") )

message("done")