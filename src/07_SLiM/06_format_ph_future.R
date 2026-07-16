# Nucella ABC analysis
## Part 3 - ABC estimation

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'foreach', 'lme4', 'abc', 'RColorBrewer'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(foreach)
library(lme4)
library(abc)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
real_All <- get(load("data/processed/SLiM/ph_ABC/real_data.Rdata"))
sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data.Rdata"))
#sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_linear.Rdata"))
