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
library(minpack.lm)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/Mcali")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# ================================================================================== #

# Load data
afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T)
# Extract just sites and env var
Mcali <- afs.Mcali[,c(7,13)] %>% distinct()

# ================================================================================== #

# Make Site factor
afs.Mcali$Site <- factor(afs.Mcali$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Subset data
afs.Mcali_sub <- afs.Mcali[,c(10,13)]

# Fit using self-starting parameters
mod <- nls(AF ~ SSlogis(mean_integrated_thk, Asym, xmid, scal), data = afs.Mcali_sub)
mod_fit <- coef(mod)

# Graph data with sigmoid curve
pdf("output/figures/SLiM/Mcali/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
ggplot(afs.Mcali, aes(x = mean_integrated_thk, y = AF, fill = Site)) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = mycolors) +
stat_function(fun = SSlogis, args = list(Asym = mod_fit["Asym"], xmid = mod_fit["xmid"], scal = mod_fit["scal"])) + theme_linedraw()
dev.off()

# ================================================================================== #


# Function

# Logistic sigmoid scaled to [-1, 1]
# where:
# z = switching point
# k = width of transition


# Selection
sel <- function(x, z, k) {
  1 / (1 + exp((x - z)/k)) - 0.5
}

# Graph
pdf("output/figures/SLiM/Mcali/Mcali_selection.pdf", width = 5, height = 5)
plot(afs.Mcali$mean_integrated_thk, sel(afs.Mcali$mean_integrated_thk, 1.92, 0.32))
dev.off()
pdf("output/figures/SLiM/Mcali/Mcali_selection_zoomed_out.pdf", width = 5, height = 5)
plot(seq(0, 5, by = 0.01), sel(seq(0, 5, by = 0.01), 1.98, 0.12))
dev.off()

# Fitness curves
fit <- function(x, z, k) {
  1 + (1 / (1 + exp((x - z)/k)) - 0.5)
}

# Graph
pdf("output/figures/SLiM/Mcali/Mcali_fitness.pdf", width = 5, height = 5)
plot(afs.Mcali$mean_integrated_thk, fit(afs.Mcali$mean_integrated_thk, 1.98, 0.12))
dev.off()

# ================================================================================== #

# Make guide file with broader range of variables

# Range of values for each parameter
thresh <- seq(1.0, 2.5, by = 0.1) #16
k <- seq(0.05, 0.20, by = 0.01) #16
mag <- c(1) #1
m <- c(0.005) #1
N <- c(5000) #1

# Make every combination of variables - 256 combos
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #

# Make guide file with broader range of variables - v2

# Range of values for each parameter
thresh <- seq(1.6, 2.2, by = 0.025) #25
k <- seq(0.1, 0.7, by = 0.05) #13
mag <- c(1) #1
m <- c(0.005) #1
N <- c(5000) #1

# Make every combination of variables - 325 combos
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file_v2.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #

# Make guide file with broader range of variables - v3

# Range of values for each parameter
thresh <- seq(1.7, 2.0, by = 0.01) #31
k <- seq(0.62, 1, by = 0.02) #
mag <- c(1) #1
m <- c(0.005) #1
N <- c(5000) #1

# Make every combination of variables - 806
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file_v3.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


# ================================================================================== #

# Make guide file with broader range of variables - v4

# Range of values for each parameter
thresh <- seq(1.7, 1.95, by = 0.01) #26
k <- seq(0.5, 1.6, by = 0.1) #11
mag <- c(1) #1
m <- c(0.0001, 0.0005, 0.0010) #3
N <- c(5000) #1

# Make every combination of variables - 936
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file_v4.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #

# Make guide file with broader range of variables - v5

# Range of values for each parameter
thresh <- seq(1.73, 2.05, by = 0.01) #
k <- seq(1, 10, by = 1) #20
mag <- c(1) #1
m <- c(0.0001, 0.0005, 0.0010) #3
N <- c(5000) #1

# Make every combination of variables - 990
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file_v5.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #

# Make guide file with broader range of variables - v6

# Range of values for each parameter
#thresh <- seq(1.85, 2.02, by = 0.01) #21
thresh <- seq(1.8, 1.84, by = 0.01) #21
k <- seq(2, 50, by = 2) #25
mag <- c(1) #1
m <- c(0.0001, 0.0010) #2
N <- c(5000) #1

# Make every combination of variables - 900 combos (250 pt2)
guide_file <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file, file = "guide_files/slim_Mcali_guide_file_v6_expand.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)
