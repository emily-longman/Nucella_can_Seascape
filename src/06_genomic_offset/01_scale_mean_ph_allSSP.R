# Scale mean pH

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load data

# Read current ph cov file
ph.cov.file <- read.table("guide_files/Baypass_ph_mean.txt")

# Load ph var in future (2090 ssp 126)
ph.cov.file.future.ssp119 <- read.table("guide_files/Baypass_ph_mean_ssp119_future.txt")
ph.cov.file.future.ssp126 <- read.table("guide_files/Baypass_ph_mean_ssp126_future.txt")
ph.cov.file.future.ssp245 <- read.table("guide_files/Baypass_ph_mean_ssp245_future.txt")
ph.cov.file.future.ssp370 <- read.table("guide_files/Baypass_ph_mean_ssp370_future.txt")
ph.cov.file.future.ssp585 <- read.table("guide_files/Baypass_ph_mean_ssp585_future.txt")

# ================================================================================== #

# Redo analyses and scale covariable

# Calc mean
ph.mean <- mean(t(ph.cov.file))
# Calc sd
ph.sd <- sd(t(ph.cov.file))

# Scale
ph.cov.file.scaled <- (ph.cov.file - ph.mean)/ph.sd

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.ssp119.scaled <- (ph.cov.file.future.ssp119 - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.future.ssp119.scaled, "guide_files/Baypass_ph_mean_ssp119_future_scaled.txt", col.names=F, row.names=F)

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.ssp126.scaled <- (ph.cov.file.future.ssp126 - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.future.ssp126.scaled, "guide_files/Baypass_ph_mean_ssp126_future_scaled.txt", col.names=F, row.names=F)

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.ssp245.scaled <- (ph.cov.file.future.ssp245 - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.future.ssp245.scaled, "guide_files/Baypass_ph_mean_ssp245_future_scaled.txt", col.names=F, row.names=F)

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.ssp370.scaled <- (ph.cov.file.future.ssp370 - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.future.ssp370.scaled, "guide_files/Baypass_ph_mean_ssp370_future_scaled.txt", col.names=F, row.names=F)

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.ssp585.scaled <- (ph.cov.file.future.ssp585 - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.future.ssp585.scaled, "guide_files/Baypass_ph_mean_ssp585_future_scaled.txt", col.names=F, row.names=F)
