# Format Bio-Oracle data (https://www.bio-oracle.org/index.php)

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'psych'))
library(data.table)
library(tidyverse)

# ================================================================================== #
# ================================================================================== #

# Load data

# Read in Bio-oracle present data
ph <- read.csv("data/raw/Bio-oracle/present/ph_baseline_2000_2018_depthsurf_taller_area.csv", header=T)
ph_future <- read.csv("data/raw/Bio-oracle/future/ph_ssp585_2020_2100_depthsurf_taller_area.csv", header=T)

# ================================================================================== #

# Remove first row (i.e. units)
ph <- ph[-1,]

# Change latitude and longitude to numeric
ph <- ph %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# ================================================================================== #

# Write table
write.csv(ph, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_taller_area.csv", row.names=FALSE)

# ================================================================================== #

# Remove first row (i.e. units)
ph_future <- ph_future[-1,]

# Change latitude and longitude to numeric
ph_future <- ph_future %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# ================================================================================== #

# Write table
write.csv(ph_future, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_ph_future_taller_area.csv", row.names=FALSE)
