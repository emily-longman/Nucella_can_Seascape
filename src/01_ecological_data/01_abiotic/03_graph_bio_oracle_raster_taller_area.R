# Graph Bio-Oracle data (https://www.bio-oracle.org/index.php)
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

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

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'terra', 'raster'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(terra)
library(raster)

# ================================================================================== #

# Figure directory
out_dir_fig <- paste("output/figures/enviro_data/Bio-oracle")
if (!dir.exists(out_dir_fig)) {dir.create(out_dir_fig)}

# ================================================================================== #

# Read in Bio-oracle present data
bio_oracle <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_taller_area.csv", header=T)

# Extract the most recent data (i.e., 2010-01-01T00:00:00Z, which represents 2010-~2020)
bio_oracle_2010 <- bio_oracle %>% 
  filter(time == "2010-01-01T00:00:00Z")

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(32, 48.5)
longitude_range <- c(-125, -118.5)

# Set study extent
study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])

# Raster resolution (focal cells of bio-oracle data are at 0.05 degree resolution)
raster_resolution <- 0.05

# Create raster layer object - specify study extent, resolution and coordinate reference system
study_raster <- raster(study_extent, res=raster_resolution, crs="+proj=longlat +datum=WGS84")

# ================================================================================== #

# Create empty raster template

# Fill dataset with NA
values(study_raster) <- NA

# ================================================================================== #

# Test graphing just one environmental variable (e.g., mean temperature)

# Set coordinates of bio-oracle data
coordinates <- cbind(bio_oracle_2010$longitude, bio_oracle_2010$latitude)

# Extract one variable - ph mean 
bio_oracle_2010_ph_mean <- bio_oracle_2010 %>% 
  dplyr::select(longitude, latitude, ph_mean)


# Rasterize ph data
ph_raster <- rasterize(coordinates, study_raster, bio_oracle_2010_ph_mean$ph_mean, fun = mean, na.rm = TRUE)

# Graph with ggplot
# Change to data frame
raster_df_ph <- as.data.frame(ph_raster, xy = TRUE, na.rm = TRUE)

# Graph pH (note: reversed colors since low pH is the stressor)
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_ph_mean_ggplot_taller.pdf",  width = 6, height = 15) 
ggplot(raster_df_ph, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-125, -122, -119)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="mean pH", breaks = c(7.92, 7.98, 8.04)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  labs(x = "Longitude", y = "Latitude") + 
  theme_linedraw(base_size = 32) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  theme(legend.title = element_text(size = 26), legend.text = element_text(size = 20), legend.position = c(0.75, 0.85), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
  #theme(plot.title = element_text(hjust=0.5))
dev.off()
