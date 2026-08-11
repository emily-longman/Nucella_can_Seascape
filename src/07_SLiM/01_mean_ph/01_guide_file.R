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
out_fig_dir <- paste("output/figures/SLiM/ph")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
afs.ph <- read.csv("data/processed/SLiM/afs.ph.g27343.BF.POD.csv", header=T)
# Extract just sites and env var
ph <- afs.ph[, c(2,8)] %>% distinct()

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
#s_v2_alt <- function(x, z, k) {
#  -1 * (2 / (1 + exp((x - z)/k)) - 1)
#}

# Graph
pdf("output/figures/SLiM/ph/pH_dist_sel_v2.pdf", width = 5, height = 5)
plot(seq(7, 9, by = 0.01), s_v2(seq(7, 9, by = 0.01), 7.99, 0.1))
dev.off()
#pdf("output/figures/SLiM/mean_ph/pH_dist_sel_v2_alt.pdf", width = 5, height = 5)
#plot(ph$ph_mean, s_v2_alt(ph$ph_mean, 7.975, 0.1))
#dev.off()


# Fitness curves
f_v2 <- function(x, z, k) {
  1 + (2 / (1 + exp((x - z)/k)) - 1)
}
#f_v2_alt <- function(x, z, k) {
#  1 + (-1 * (2 / (1 + exp((x - z)/k)) - 1))
#}

# Graph
pdf("output/figures/SLiM/ph/pH_dist_v2.pdf", width = 5, height = 5)
plot(ph$ph_mean, f_v2(ph$ph_mean, 7.975, 0.01))
dev.off()
#pdf("output/figures/SLiM/mean_ph/pH_dist_v2_alt.pdf", width = 5, height = 5)
#plot(ph$ph_mean, f_v2_alt(ph$ph_mean, 7.975, 0.01))
#dev.off()

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
pdf("output/figures/SLiM/ph/pH_sel_magnitude.pdf", width = 5, height = 5)
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
pdf("output/figures/SLiM/ph/pH_dist_sel_magnitude.pdf", width = 5, height = 5)
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


####

# More vars 3

# Range of values for each parameter
thresh <- seq(7.94, 8.02, by=0.005) #17
k <- seq(0.05, 0.25, by = 0.01) #5
#mag <- c(1) #1 - pt1 
mag <- c(2) #1
m <- c(0.001) #1
N <- c(2500) #1

# Make every combination of variables - 90 combos
guide_file_morevars3_pt2 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars3_pt2, file = "guide_files/slim_ph_guide_file_morevars3_pt2.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


###

# More vars 4

# Range of values for each parameter
thresh <- seq(7.98, 8.01, by=0.001) #31
#thresh <- seq(7.98, 8.01, by=0.0025) #13
#k <- seq(0.05, 0.1, by = 0.01) #6
k <- seq(0.11, 0.15, by = 0.01) #6
mag <- c(1) #1
#m <- seq(0.001, 0.005, by = 0.001) #5
m <- seq(0.006, 0.010, by = 0.001) #5
N <- c(5000) #1

# Make every combination of variables - 775 #930 combos
guide_file_morevars4 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars4, file = "guide_files/slim_ph_guide_file_morevars4_biggerrange.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


###

# More vars 5

# Range of values for each parameter
thresh <- seq(7.98, 7.995, by=0.001) #16
k <- seq(0.04, 0.28, by = 0.04) #11
mag <- c(1) #1
m <- seq(0.005, 0.040, by = 0.005) #8
N <- c(5000) #1

# Make every combination of variables - 896 combos
guide_file_morevars5 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars5, file = "guide_files/slim_ph_guide_file_morevars5.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


s <- function(x, z, k, mag) {
  mag / (1 + exp((x - z)/k)) - (mag/2)
}

# Graph
pdf("output/figures/SLiM/ph_ABC/pH_sel_magnitude.pdf", width = 5, height = 5)
plot(seq(7.8, 8.1, by = 0.01), s(seq(7.8, 8.1, by = 0.01), 7.987, 0.15, 1))
dev.off()


###

# More vars 6

# Range of values for each parameter
thresh <- seq(7.98, 7.99, by = 0.001) #11
k <- seq(0.01, 0.07, by = 0.005) #13
mag <- c(1) #1
m <- seq(0.001, 0.005, by = 0.001) #5
N <- c(5000) #1

# Make every combination of variables - 715 combos
guide_file_morevars6 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars6, file = "guide_files/slim_ph_guide_file_morevars6.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

###

# More vars 7

# Range of values for each parameter
thresh <- seq(7.980, 7.992, by = 0.0005) #11
#k <- seq(0.01, 0.10, by = 0.005) #13
k <- seq(0.105, 0.15, by = 0.005) #13
mag <- c(1) #1
m <- c(0.001) #1
N <- c(5000) #1

# Make every combination of variables - 475 combos (250 combos pt2)
guide_file_morevars7 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars7, file = "guide_files/slim_ph_guide_file_morevars7_expand.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

###

# More vars 8

# Range of values for each parameter
thresh <- seq(7.980, 7.992, by = 0.0005) #11
k <- seq(0.01, 0.10, by = 0.005) #13
mag <- c(1) #1
m <- c(0.00001) #1
N <- c(5000) #1

# Make every combination of variables - 475 combos
guide_file_morevars8 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars8, file = "guide_files/slim_ph_guide_file_morevars8.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)


###

# More vars 9

# Range of values for each parameter
thresh <- seq(7.980, 7.992, by = 0.0005) #11
k <- seq(0.01, 0.10, by = 0.005) #13
mag <- c(1) #1
m <- c(0.0001) #1
N <- c(5000) #1

# Make every combination of variables - 475 combos
guide_file_morevars9 <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars9, file = "guide_files/slim_ph_guide_file_morevars9.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# More vars 9 expand

# Range of values for each parameter
thresh <- seq(7.980, 7.992, by = 0.0005) #11
#k <- seq(0.01, 0.10, by = 0.005) #13
k <- seq(0.11, 0.20, by = 0.005) #13
mag <- c(1) #1
m <- c(0.0001) #1
N <- c(5000) #1

# Make every combination of variables - 475 combos
guide_file_morevars9_expand <- expand.grid(thresh, k, mag, m, N)

# Write table
write.table(guide_file_morevars9_expand, file = "guide_files/slim_ph_guide_file_morevars9_expand.txt", sep = "\t", quote = FALSE, row.names=F, col.names=F)

# ================================================================================== #
# ================================================================================== #


# Graph and cal sigmoid for real data

# Make Site factor
topsnp$Site <- factor(topsnp$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Graph real
pdf("output/figures/SLiM/mean_ph/real_AFs_topsnp_sigmoid.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = AF, y = ph_mean, fill = Site)) + geom_point(size = 3, shape = 21) +  scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()
pdf("output/figures/SLiM/mean_ph/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = ph_mean, y = AF, fill = Site)) + geom_point(size = 3, shape = 21) +  scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()

# Subset data
topsnp_sub <- topsnp[,c(5,8)]

# Fit using self-starting parameters
mod <- nls(AF ~ SSlogis(ph_mean, Asym, xmid, scal), data = topsnp_sub)
mod_fit <- coef(mod)

pdf("output/figures/SLiM/mean_ph/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
plot(topsnp_sub$ph_mean, topsnp_sub$AF, pch = 20)
curve(SSlogis(x, mod_fit["Asym"], mod_fit["xmid"], mod_fit["scal"]), lwd = 2, col = 'lightblue', add = TRUE)
dev.off()

pdf("output/figures/SLiM/mean_ph/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = ph_mean, y = AF, fill = Site)) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = mycolors) +
stat_function(fun = SSlogis, args = list(Asym = mod_fit["Asym"], xmid = mod_fit["xmid"], scal = mod_fit["scal"])) + theme_linedraw()
dev.off()
