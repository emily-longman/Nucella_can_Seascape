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

# Present day data

# Read in Bio-oracle data
sst <- read.csv("data/raw/Bio-oracle/present/thetao_baseline_2000_2019_depthsurf_larger_area.csv", header=T)
chl <- read.csv("data/raw/Bio-oracle/present/chl_baseline_2000_2018_depthsurf_larger_area.csv", header=T)
o2 <- read.csv("data/raw/Bio-oracle/present/o2_baseline_2000_2018_depthsurf_larger_area.csv", header=T)
ph <- read.csv("data/raw/Bio-oracle/present/ph_baseline_2000_2018_depthsurf_larger_area.csv", header=T)
so <- read.csv("data/raw/Bio-oracle/present/so_baseline_2000_2019_depthsurf_larger_area.csv", header=T)

# Read in future Bio-oracle data
sst_future <- read.csv("data/raw/Bio-oracle/future/thetao_ssp585_2020_2100_depthsurf_larger_area.csv", header=T)
chl_future <- read.csv("data/raw/Bio-oracle/future/chl_ssp585_2020_2100_depthsurf_larger_area.csv", header=T)
o2_future <- read.csv("data/raw/Bio-oracle/future/o2_ssp585_2020_2100_depthsurf_larger_area.csv", header=T)
ph_future <- read.csv("data/raw/Bio-oracle/future/ph_ssp585_2020_2100_depthsurf_larger_area.csv", header=T)
so_future <- read.csv("data/raw/Bio-oracle/future/so_ssp585_2020_2100_depthsurf_larger_area.csv", header=T)

# ================================================================================== #

# Combine datasets into one large dataset
bio_oracle <- reduce(list(sst, chl, o2, ph, so), full_join, by = c("time", "latitude", "longitude"))

# Remove first row (i.e. units)
bio_oracle <- bio_oracle[-1,]

# Change latitude and longitude to numeric
bio_oracle <- bio_oracle %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# ================================================================================== #

# Write table
write.csv(bio_oracle, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_larger_area.csv", row.names=FALSE)

# ================================================================================== #
# ================================================================================== #

# Combine datasets into one large dataset
bio_oracle_future <- reduce(list(sst_future, chl_future, o2_future, ph_future, so_future), full_join, by = c("time", "latitude", "longitude"))

# Remove first row (i.e. units)
bio_oracle_future <- bio_oracle_future[-1,]

# Change latitude and longitude to numeric
bio_oracle_future <- bio_oracle_future %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# ================================================================================== #

# Write table
write.csv(bio_oracle_future, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_future_larger_area.csv", row.names=FALSE)
