# Format pH future Bio-oracle data

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
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]

# ================================================================================== #

# Generate output directories

# Bio-Oracle directory
out_dir_future <- paste("data/raw/Bio-oracle/future")
if (!dir.exists(out_dir_future)) {dir.create(out_dir_future)}

# Data SliM directory
out_data_dir <- paste("data/processed/SLiM/ph_future")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# Extract future day data 2020-2100

# Define time, lat, and long (set lat and long to fully encompass all sites)
time = c('2020-01-01T00:00:00Z', '2090-01-01T00:00:00Z')
latitude = c(34, 45)
longitude = c(-120, -125)

# Set constraints
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

# Set dataset IDs
pH_SSP585 <- "ph_ssp585_2020_2100_depthsurf"

# Specify datasets and summary statistics to download 
datasets <- list(list(dataset_id = pH_SSP585, variables = c("ph_mean"),
       constraints = constraints, fmt = "csv", directory = out_dir_future))

# Download dataset
for (dataset in datasets) {
  
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  # List files before downloading
  files_before <- list.files(out_dir_future, full.names = TRUE)
  # Download dataset
  download_layers(dataset_id, variables = variables, constraints = constraints, fmt = "csv", directory = out_dir_future)
  # List files after downloading
  files_after <- list.files(out_dir_future, full.names = TRUE)
  # Identify the newly downloaded file
  new_file <- setdiff(files_after, files_before)
  
  # Rename the file to match the dataset ID
  if (length(new_file) == 1) {
    new_name <- file.path(out_dir_future, paste0(dataset_id, "_decadal.csv"))
    file.rename(new_file, new_name)
    message("Rename ", new_file, " to ", new_name)
  } else {
    message("No new file found for ", dataset_id, " or multiple new files detected.")
  }
}

# ================================================================================== #

# Read in Bio-oracle data
ph_ssp585 <- read.csv("data/raw/Bio-oracle/future/ph_ssp585_2020_2100_depthsurf_decadal.csv", header=T)

# ================================================================================== #

# Extract data for just sites

# Remove first row (i.e. units)
ph_ssp585 <- ph_ssp585[-1,]

# Change latitude and longitude to numeric
ph_ssp585 <- ph_ssp585 %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# Filter bio-oracle data for the 19 sites (will set NA for any lat long combo not specified, then will filter out those rows)
# Note: latitude and longitude need to be rounded to nearest 0.025 or 0.075
# Note: during rounding make sure you don't set the location as a spot on land
ph_ssp585_sites <- ph_ssp585 %>%
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

# ================================================================================== #

# Add decade term
ph_ssp585_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_lm <- foreach(i=unique(ph_ssp585_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp585_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp585_sites$location <- factor(ph_ssp585_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_lm$Site <- factor(ph_future_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_lm <- ph_future_lm[order(ph_future_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_lm, "data/processed/SLiM/ph_future/ph_future_lm.csv", row.names=F)

# ================================================================================== #

# Extract last decade 
ph_ssp585_sites_2090_tmp <- ph_ssp585_sites[which(ph_ssp585_sites$decade == 75),]
# Reverse order
ph_ssp585_sites_2090 <- ph_ssp585_sites_2090_tmp[nrow(ph_ssp585_sites_2090_tmp):1,]

write.csv(ph_ssp585_sites_2090, "data/processed/SLiM/ph_future/ph_future_2090_sites.csv", row.names=F)

# ================================================================================== #

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle.pdf", width = 6.5, height = 6)
ggplot(ph_ssp585_sites, aes(x = decade, y = ph_mean, color = location)) + geom_point(size = 3, shape = 16, alpha = 0.9) +
    geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=7.987, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()


# For visualization purposes make x-axis go from 2020-2100

# Calc pred pH adjusted for 2020-2100
pred <- ph_future_lm %>%rowwise() %>%
do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Threshold from simulations
thresh = 7.986

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_alt.pdf", width = 6.5, height = 6)
ggplot(ph_ssp585_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 3, shape = 16, alpha = 0.9) +
    geom_line(data = pred, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# 4 population subset

# Graph subset of poplations
ph_ssp585_sites_sub <- ph_ssp585_sites %>% filter(location == "FC" | location == "CBL" | location == "BMR" | location == "STR")
ph_ssp585_sites_sub$location <- factor(ph_ssp585_sites_sub$location, levels = c("FC", "CBL", "BMR", "STR"))
pred_sub <- pred %>% filter(Site == "FC" | Site == "CBL" | Site == "BMR" | Site == "STR")
pred_sub$Site <- factor(pred_sub$Site, levels = c("FC", "CBL", "BMR", "STR"))

# Graph subset
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_subset.pdf", width = 7, height = 5.5)
ggplot(ph_ssp585_sites_sub, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 3, shape = 16, alpha = 0.9) +
    geom_line(data = pred_sub, aes(x = year, y = ph, color = Site), linewidth = 2.5)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors_sub) + scale_color_manual(values = mycolors_sub) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()


# ================================================================================== #

# When will FC be below the threshold

# Subset
FC <- pred[which(pred$Site=="FC"),]

# Calc year when FC will be below thresh
min(FC$year[which(FC$ph<thresh)])
