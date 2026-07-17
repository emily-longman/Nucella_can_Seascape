# Create guide file of parameters

# Clear memory
rm(list=ls())

# Stop exponential
options(scipen = 999)

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
plot(seq(7, 9, by = 0.01), s_v2(seq(7, 9, by = 0.01), 7.99, 0.1))
dev.off()
pdf("output/figures/SLiM/pH_dist_sel_v2_alt.pdf", width = 5, height = 5)
plot(ph$ph_mean, s_v2_alt(ph$ph_mean, 7.975, 0.1))
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

# Make guide file with broader range of variables

# Range of values for each parameter
thresh <- seq(7.7, 8.4, by=0.05) #19
k_1 <- c(0, 0.001, 0.01, 0.1) #4
k_2 <- c(0, 0.001, 0.01, 0.1) #4
m <- c(0.01, 0.001, 0.0001, 0.00001) #4

# Make every combination of variables - 960 combos
guide_file_expandedparam <- expand.grid(thresh, k_1, k_2, m)

# Write table
write.table(guide_file_expandedparam, file = "guide_files/slim_ph_guide_file_expandedparam.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


# ================================================================================== #
# ================================================================================== #


# Other ideas
# Selection
s <- function(x, z, k, mag) {
  mag / (1 + exp((x - z)/k)) - (mag/2)
}

# Graph
pdf("output/figures/SLiM/pH_dist_sel_magnitude.pdf", width = 5, height = 5)
plot(ph$ph_mean, s(ph$ph_mean, 7.975, 0.01, 4))
dev.off()

# Make guide file with broader range of variables

# Range of values for each parameter
thresh <- seq(7.94, 8.02, by=0.01) #9
k <- c(0.12, 0.10, 0.08, 0.06) #4
mag <- c(1, 2, 3) #3 -- note: 2 matches what my original sel function had
m <- c(0.01, 0.001, 0.0001) #3
N <- c(1000, 2500, 5000) #3

# Make every combination of variables - 972 combos
guide_file_morevars <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars, file = "guide_files/slim_ph_guide_file_morevars.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #
# ================================================================================== #

# More vars 2

# Selection
s <- function(x, z, k, mag) {
  mag / (1 + exp((x - z)/k)) - (mag/2)
}

# Graph
pdf("output/figures/SLiM/pH_dist_sel_magnitude.pdf", width = 5, height = 5)
plot(ph$ph_mean, s(ph$ph_mean, 7.99, 0.06, 1.5))
dev.off()

# Make guide file with broader range of variables

# Range of values for each parameter
thresh <- seq(7.98, 8.005, by=0.005) #6
k <- c(0.08, 0.07, 0.06, 0.05, 0.04) #5
mag <- c(1, 1.5, 2) #3
m <- c(0.01, 0.005,  0.001) #3
N <- c(4000, 5000, 6000) #3

# Make every combination of variables - 810 combos
guide_file_morevars2 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars2, file = "guide_files/slim_ph_guide_file_morevars2.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


###

# More vars 3

# Range of values for each parameter
thresh <- seq(7.94, 8.02, by=0.005) #17
k <- seq(0.05, 0.25, by = 0.01) #5
#mag <- c(1) #1
mag <- c(2) #1
m <- c(0.001) #1
N <- c(2500) #1

# Make every combination of variables - 90 combos
guide_file_morevars3_pt2 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars3_pt2, file = "guide_files/slim_ph_guide_file_morevars3_pt2.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)
