# Create guide file of parameters

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

# OLD
#f <- function(x, z, k) {
#  1 + (2 / (1 + exp(-(x - z)/k)) - 1)
#}
#f_alt <- function(x, z, k) {
#  1 + (-1 * (2 / (1 + exp(-(x - z)/k)) - 1))
#}
# OLD
#pdf("output/figures/SLiM/pH_dist.pdf", width = 5, height = 5)
#plot(ph$ph_mean, f(ph$ph_mean, 7.975, 0.001))
#dev.off()

# Selection
s_v2 <- function(x, z, k) {
  2 / (1 + exp((x - z)/k)) - 1
}
s_v2_alt <- function(x, z, k) {
  -1 * (2 / (1 + exp((x - z)/k)) - 1)
}

# Graph
pdf("output/figures/SLiM/pH_dist_sel_v2.pdf", width = 5, height = 5)
plot(ph$ph_mean, s_v2(ph$ph_mean, 7.975, 0.0001))
dev.off()
pdf("output/figures/SLiM/pH_dist_sel_v2_alt.pdf", width = 5, height = 5)
plot(ph$ph_mean, s_v2_alt(ph$ph_mean, 7.975, 0.01))
dev.off()


# Fitness curves
f_v2 <- function(x, z, k) {
  1 + (2 / (1 + exp((x - z)/k)) - 1)
}
f_v2_alt <- function(x, z, k) {
  1 + (-1 * (2 / (1 + exp((x - z)/k)) - 1))
}

# Graph
pdf("output/figures/SLiM/pH_dist_v2.pdf", width = 5, height = 5)
plot(ph$ph_mean, f_v2(ph$ph_mean, 7.975, 0.01))
dev.off()
pdf("output/figures/SLiM/pH_dist_v2_alt.pdf", width = 5, height = 5)
plot(ph$ph_mean, f_v2_alt(ph$ph_mean, 7.975, 0.01))
dev.off()

# Good options for k: 0.1, 0.01, 0.001

# ================================================================================== #

# Create guide file
ph_min <- plyr::round_any(min(ph$ph_mean), 0.01, f=ceiling)
ph_max <- plyr::round_any(max(ph$ph_mean), 0.01, f=floor)
k <- c(0.001, 0.01, 0.1)

# Cross-join function in data table
guide_file_tmp <- as.data.frame(CJ(seq(ph_min, ph_max, by=0.01), k))
# Rename for each k
guide_file_tmp1 <- guide_file_tmp %>% rename(thresh = V1, k_1 = k)
guide_file_tmp2 <- guide_file_tmp %>% rename(thresh = V1, k_2 = k)

# Make first 1/3 of table
guide_file_pt1 <- data.table(seq(ph_min, ph_max, by=0.01), 0, 0) %>% rename(thresh = V1, k_1 = V2, k_2 = V3)
# Make second 1/3 of table
guide_file_pt2 <- rbind(data.table(guide_file_tmp1, 0),
                    data.table(guide_file_tmp1, rep(k, 10))) %>% rename(k_2 = V2)
# Make third 1/3 of table
guide_file_pt3 <- guide_file_tmp2 %>% mutate(k_1 = 0, .after = thresh)

# Bind
guide_file <- rbind(guide_file_pt1, guide_file_pt2, guide_file_pt3)

# ================================================================================== #

# Write table
write.table(guide_file, file = "guide_files/slim_ph_guide_file.txt", sep = "\t", row.names=F, col.names=F)


# ================================================================================== #
# ================================================================================== #

# Make guide file with broader range of variables

# Range of values for each parameter
thresh <- seq(7.7, 8.4, by=0.05) #19
k_1 <- c(0, 0.001, 0.01, 0.1) #4
k_2 <- c(0, 0.001, 0.01, 0.1) #4
m <- c(0.01, 0.001, 0.0001, 0.0001) #4

# Make every combination of variables - 960 combos
guide_file_expandedparam <- expand.grid(thresh, k_1, k_2, m)

# Write table
write.table(guide_file_expandedparam, file = "guide_files/slim_ph_guide_file_expandedparam.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

