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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'psych', 'viridis'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# ================================================================================== #

# Set colors
viridiscolors <- viridis(n=19)

# ================================================================================== #

# Present day data

# Read in Bio-oracle data
ph <- read.csv("data/raw/Bio-oracle/present/pH/ph_baseline_2000_2018_depthsurf.csv", header=T)

# ================================================================================== #

# Remove first row (i.e. units)
ph <- ph[-1,]

# Change latitude and longitude to numeric
ph <- ph %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude))

# Filter bio-oracle data for the 19 sites (will set NA for any lat long combo not specified, then will filter out those rows)
# Note: latitude and longitude need to be rounded to nearest 0.025 or 0.075
# Note: during rounding make sure you don't set the location as a spot on land
ph_sites <- ph %>%
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
str(ph_sites)

# ================================================================================== #

# Extract data for most recent decade
ph_sites_2010 <- ph_sites %>% 
filter(time == "2010-01-01T00:00:00Z")

# ================================================================================== #

# Set location as factor
ph_sites_2010$location <- factor(ph_sites_2010$location, 
    levels = rev(c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")))

# Add column that highlights N, and S
ph_sites_2010 <- ph_sites_2010 %>% mutate(shape = case_when(location %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   location %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))

# Graph
pdf("output/figures/enviro/Bio-oracle/pH_sites.pdf", width = 6, height = 7.18)
ggplot(ph_sites_2010, aes(x = ph_mean, y = location, shape = shape)) +
  geom_linerange(data = ph_sites_2010, aes(xmin = ph_min, xmax = ph_max)) +
  scale_shape_manual(values = c(21, 23)) +
  geom_point(size=6, color = "black", fill = rev(viridiscolors)) +
  labs(x = "pH", y = "") + xlim(7.7, 8.2) +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()