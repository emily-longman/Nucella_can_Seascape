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

# ================================================================================== #

# Format data

# Create SNP_id column
glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Create windows

# Define window and step size
win.bp <- 1e5
step.bp <- 5e8

# How many SNPs are on each contig:
SNPS_density <- glm.model.collated %>% group_by(chr) %>% summarize(n=n())
# Graph
pdf("output/figures/GEA/glms/glms_window_summary/glm_pval_density.pdf", width = 8, height = 8)
ggplot(SNPS_density, aes(x=n))+ geom_density()
dev.off()
# Use this information to determine level to filter for number of SNPs in a given window

# Generate windows (note: only windows with the number of SNPs in that window >= 5)
wins <- foreach(chr.i=unique(glm.model.collated$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  
                  tmp <- glm.model.collated %>%
                    filter(chr == chr.i)
                  
                  nSNPs=dim(tmp)[1]
                  
                  if(nSNPs >= 5){
                    o =
                      data.table(chr=chr.i,
                                 nSNPs=dim(tmp)[1],
                                 start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                                 end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                    return(o)
                    
                  }   
                  else {message("fails nSNPs filter")}
                }

# Add window index
wins[,i:=1:dim(wins)[1]]

# Check dimensions
dim(wins)

# Save windows
save(wins, file="data/processed/GEA/glms/glms_window_summary/windows.RData")

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


# Window summarization for each  environmental var
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