# Format Mcali data

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'devtools', 'magrittr, 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(tidyverse)
library(devtools)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Data directory
out_data_dir <- paste("data/processed/SLiM/Mcali_future")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# 
