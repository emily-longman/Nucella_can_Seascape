# Filter GLM so only real

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
#install.packages(c('tidyverse', 'dplyr'))
library(tidyverse)
library(dplyr)

# ================================================================================== #

# Load GLM data (glm.model.collated.filt)
load("data/processed/GEA/glms/glms_per_var/glm.collated_ph_mean_filt.Rdata")

# ================================================================================== #

# Filter
glm.model.collated.filt.real <- glm.model.collated.filt %>% filter(data=="real")

# ================================================================================== #

# Save
save(glm.model.collated.filt.real, file="data/processed/GEA/glms/glms_per_var/glm.collated_ph_mean_filt_real.Rdata")