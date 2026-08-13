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

# Load ph var in future (2090 ssp 585)
ph.cov.file.future <- read.table("guide_files/Baypass_ph_mean_future.txt")

# ================================================================================== #

# Redo analyses and scale covariable

# Calc mean
ph.mean <- mean(t(ph.cov.file))
# Calc sd
ph.sd <- sd(t(ph.cov.file))

# Scale
ph.cov.file.scaled <- (ph.cov.file - ph.mean)/ph.sd

# Write file
write.table(ph.cov.file.scaled, "guide_files/Baypass_ph_mean_scaled.txt", col.names=F, row.names=F)

# Scale future with respect to the mean and variance of the original covariable values
ph.cov.file.future.scaled <- (ph.cov.file.future - ph.mean)/ph.sd

# Write file
write.table(ph.cov.file.future.scaled, "guide_files/Baypass_ph_mean_future_scaled.txt", col.names=F, row.names=F)
