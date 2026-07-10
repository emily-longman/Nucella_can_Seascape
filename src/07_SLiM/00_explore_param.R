# Test

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
#install.packages(c('data.table', 'tidyverse', 'foreach'))
library(data.table)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/SLiM")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
afs.ph <- read.csv("data/processed/SLiM/afs.ph.g27343.BF.POD.csv", header=T)
afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T)

# Extract just sites and env var
ph <- afs.ph[, c(2,8)] %>% distinct()
Mcali <- afs.Mcali[,c(7,13)] %>% distinct()

# ================================================================================== #

# Function

# Logistic sigmoid scaled to [-1, 1]
# where:
# z = switching point
# k = width of transition
f <- function(x, z, k) {
  1 + (2 / (1 + exp(-(x - z)/k)) - 1)
}

pdf("output/figures/SLiM/pH_dist.pdf", width = 5, height = 5)
plot(ph$ph_mean, f(ph$ph_mean, 7.975, 0.1))
dev.off()

# Good options for k: 0.1, 0.01, 0.001

# ================================================================================== #

# Create guide file
ph_min <- plyr::round_any(min(ph$ph_mean), 0.01, f=ceiling)
ph_max <- plyr::round_any(max(ph$ph_mean), 0.01, f=floor)
k <- c(0.001, 0.01, 0.1)
#ks <- data.table(rep(k,2),c(rep(0, 3), k))

# Cross-join function in data table
guide_file_tmp <- as.data.frame(CJ(seq(ph_min, ph_max, by=0.01), k))
guide_file <- rbind(data.table(guide_file_tmp, 0),
                    data.table(guide_file_tmp, rep(k, 10)))

# Write table
write.table(guide_file, file = "guide_files/slim_ph_guide_file.txt", sep = "\t", row.names=F, col.names=F)




# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# TESTING

# Get test distribution of data
a = runif(18, 1.4, 2.2)

quantile(a)

selection_guess = function(s, x, thresh){
  i = (s*x)+thresh
  return(i)
}