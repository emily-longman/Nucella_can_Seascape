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
library(ggplot2)
library(RColorBrewer)
library(psych)

# ================================================================================== #
# ================================================================================== #

# Present day data

# Read in Bio-oracle data
sst <- read.csv("data/raw/Bio-oracle/present/thetao_baseline_2000_2019_depthsurf.csv", header=T)
chl <- read.csv("data/raw/Bio-oracle/present/chl_baseline_2000_2018_depthsurf.csv", header=T)
o2 <- read.csv("data/raw/Bio-oracle/present/o2_baseline_2000_2018_depthsurf.csv", header=T)
ph <- read.csv("data/raw/Bio-oracle/present/ph_baseline_2000_2018_depthsurf.csv", header=T)
so <- read.csv("data/raw/Bio-oracle/present/so_baseline_2000_2019_depthsurf.csv", header=T)

# ================================================================================== #

# Combine datasets into one large dataset
bio_oracle <- reduce(list(sst, chl, o2, ph, so), full_join, by = c("time", "latitude", "longitude"))

# Remove first row (i.e. units)
bio_oracle <- bio_oracle[-1,]

# Change latitude and longitude to numeric
bio_oracle <- bio_oracle %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# Filter bio-oracle data for the 19 sites (will set NA for any lat long combo not specified, then will filter out those rows)
# Note: latitude and longitude need to be rounded to nearest 0.025 or 0.075
# Note: during rounding make sure you don't set the location as a spot on land
bio_oracle_sites <- bio_oracle %>%
  mutate(location = case_when(
      latitude == 43.325 & longitude == -124.425 ~ "ARA",
      latitude == 38.325 & longitude == -123.075 ~ "BMR",
      latitude == 42.825 & longitude == -124.575 ~ "CBL",
      latitude == 44.825 & longitude == -124.075 ~ "FC",
      latitude == 38.525 & longitude == -123.275 ~ "FR",
      latitude == 35.275 & longitude == -120.925 ~ "HZD",
      latitude == 39.625 & longitude == -123.825 ~ "KH",
      latitude == 34.875 & longitude == -120.675 ~ "OCT",
      latitude == 35.675 & longitude == -121.325 ~ "PB",
      latitude == 37.175 & longitude == -122.375 ~ "PGP",
      latitude == 36.525 & longitude == -121.975 ~ "PL",
      latitude == 41.775 & longitude == -124.275 ~ "PSG",
      latitude == 35.725 & longitude == -121.325 ~ "PSN",
      latitude == 36.425 & longitude == -121.925 ~ "SBR",
      latitude == 44.225 & longitude == -124.125 ~ "SH",
      latitude == 44.525 & longitude == -124.075 ~ "SLR",
      latitude == 40.025 & longitude == -124.075 ~ "STC",
      latitude == 34.725 & longitude == -120.625 ~ "STR",
      latitude == 39.275 & longitude == -123.825 ~ "VD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(location))

# Check structure
str(bio_oracle_sites)

# ================================================================================== #

# Extract data for most recent decade
bio_oracle_sites_2010 <- bio_oracle_sites %>% 
filter(time == "2010-01-01T00:00:00Z")

# ================================================================================== #

# Write table
write.csv(bio_oracle, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle.csv", row.names=FALSE)
write.csv(bio_oracle_sites, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites.csv", row.names=FALSE)
write.csv(bio_oracle_sites_2010, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", row.names=FALSE)

# ================================================================================== #
# ================================================================================== #

# Future data

# Read in Bio-oracle data
sst_ssp585 <- read.csv("data/raw/Bio-oracle/future/thetao_ssp585_2020_2100_depthsurf.csv", header=T)
chl_ssp585 <- read.csv("data/raw/Bio-oracle/future/chl_ssp585_2020_2100_depthsurf.csv", header=T)
o2_ssp585 <- read.csv("data/raw/Bio-oracle/future/o2_ssp585_2020_2100_depthsurf.csv", header=T)
ph_ssp585 <- read.csv("data/raw/Bio-oracle/future/ph_ssp585_2020_2100_depthsurf.csv", header=T)
so_ssp585 <- read.csv("data/raw/Bio-oracle/future/so_ssp585_2020_2100_depthsurf.csv", header=T)

# ================================================================================== #

# Combine datasets into one large dataset
bio_oracle_ssp585 <- reduce(list(sst_ssp585, chl_ssp585, o2_ssp585, ph_ssp585, so_ssp585), full_join, by = c("time", "latitude", "longitude"))

# Remove first row (i.e. units)
bio_oracle_ssp585 <- bio_oracle_ssp585[-1,]

# Change latitude and longitude to numeric
bio_oracle_ssp585 <- bio_oracle_ssp585 %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# Filter bio-oracle data for the 19 sites (will set NA for any lat long combo not specified, then will filter out those rows)
# Note: latitude and longitude need to be rounded to nearest 0.025 or 0.075
# Note: during rounding make sure you don't set the location as a spot on land
bio_oracle_ssp585_sites <- bio_oracle_ssp585 %>%
  mutate(location = case_when(
      latitude == 43.325 & longitude == -124.425 ~ "ARA",
      latitude == 38.325 & longitude == -123.075 ~ "BMR",
      latitude == 42.825 & longitude == -124.575 ~ "CBL",
      latitude == 44.825 & longitude == -124.075 ~ "FC",
      latitude == 38.525 & longitude == -123.275 ~ "FR",
      latitude == 35.275 & longitude == -120.925 ~ "HZD",
      latitude == 39.625 & longitude == -123.825 ~ "KH",
      latitude == 34.875 & longitude == -120.675 ~ "OCT",
      latitude == 35.675 & longitude == -121.325 ~ "PB",
      latitude == 37.175 & longitude == -122.375 ~ "PGP",
      latitude == 36.525 & longitude == -121.975 ~ "PL",
      latitude == 41.775 & longitude == -124.275 ~ "PSG",
      latitude == 35.725 & longitude == -121.325 ~ "PSN",
      latitude == 36.425 & longitude == -121.925 ~ "SBR",
      latitude == 44.225 & longitude == -124.125 ~ "SH",
      latitude == 44.525 & longitude == -124.075 ~ "SLR",
      latitude == 40.025 & longitude == -124.075 ~ "STC",
      latitude == 34.725 & longitude == -120.625 ~ "STR",
      latitude == 39.275 & longitude == -123.825 ~ "VD",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(location))

# Check structure
str(bio_oracle_ssp585_sites)

# ================================================================================== #

# Extract only the last decade of data (i.e., )
bio_oracle_sites_ssp585_2090 <- bio_oracle_ssp585_sites %>% 
filter(time == "2090-01-01T00:00:00Z")

# ================================================================================== #

# Write table
write.csv(bio_oracle_ssp585, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_ssp585.csv", row.names=FALSE)
write.csv(bio_oracle_ssp585_sites, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_ssp585_sites.csv", row.names=FALSE)
write.csv(bio_oracle_sites_ssp585_2090, "data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_ssp585_2090.csv", row.names=FALSE)
