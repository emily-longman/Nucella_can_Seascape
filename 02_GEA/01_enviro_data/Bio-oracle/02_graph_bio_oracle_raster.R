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

# Read in Bio-oracle present data
bio_oracle <- read.csv("data/processed/GEA/data/bio_oracle.csv", header=T)

# Extract the most recent data (i.e., 2010-01-01T00:00:00Z, which represents 2010-2020)
bio_oracle_2010 <- bio_oracle %>% 
  filter(time == "2010-01-01T00:00:00Z")

# Read in Bio-oracle future data
bio_oracle_ssp585 <- read.csv("data/processed/GEA/data/bio_oracle_ssp585.csv", header=T)

# Extract only the last decade of data (i.e., )
bio_oracle_ssp585_2090 <- bio_oracle_ssp585 %>% 
  filter(time == "2090-01-01T00:00:00Z")

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(34, 45)
longitude_range <- c(-125, -120)

# Set study extent
study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])

# Raster resolution (focal cells of bio-oracle data are at 0.05 degree resolution)
raster_resolution <- 0.05

# Create raster layer object - specify study extent, resolution and coordinate reference system
study_raster <- raster(study_extent, res=raster_resolution, crs="+proj=longlat +datum=WGS84")

# Set color palette 
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(1000))

# ================================================================================== #

# Create empty raster template

# Fill dataset with NA
values(study_raster) <- NA

# Graph empty template
pdf("output/figures/GEA/enviro/Bio-oracle/Raster_template.pdf", width = 5, height = 5)
plot(study_raster, main = "West Coast Raster Template")
dev.off()

# ================================================================================== #

# Test graphing just one environmental variable (e.g., mean temperature)

# Set coordinates of bio-oracle data
coordinates <- cbind(bio_oracle_2010$longitude, bio_oracle_2010$latitude)

# Extract one variable - temperature 
bio_oracle_2010_temp_mean <- bio_oracle_2010 %>% 
  dplyr::select(longitude, latitude, thetao_mean)

# Rasterize temperature data
temp_raster <- rasterize(coordinates, study_raster, bio_oracle_2010_temp_mean$thetao_mean, fun = mean, na.rm = TRUE)

# Look at structure
str(temp_raster)

# Graph temperature raster
pdf("output/figures/GEA/enviro/Bio-oracle/Test_Raster_bio-oracle_temp_mean.pdf", width = 3.25, height = 5)
plot(temp_raster, col = mycolors, axes = TRUE, box = FALSE,
xlim = c(-130, 115), ylim = c(33, 46), 
xlab="Longitude", ylab="Latitude", main = "Rasterized thetao_mean")
dev.off()

# Write raster tif file
writeRaster(temp_raster, filename = "output/figures/GEA/enviro/Bio-oracle/Test_biooracle_thetao_mean_test_raster.tif", format = "GTiff")

# Graph with ggplot
# Change to data frame
raster_df <- as.data.frame(temp_raster, xy = TRUE, na.rm = TRUE)
# Graph 
pdf("output/figures/GEA/enviro/Bio-oracle/Test_Raster_bio-oracle_temp_mean_ggplot.pdf", width = 3.5, height = 5)
ggplot(raster_df, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_fill_gradientn(colours=brewer.pal(5, "RdBu")) +
  #scale_fill_viridis_c() +  # Color scale for the raster values
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  ggtitle("Rasterized thetao_mean")  +
  theme_void() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Create rasters for present day data

# Environmental variables
env_variables <- c("thetao_max",   "thetao_min",   "thetao_range", "thetao_mean", "chl_mean", "o2_mean", "ph_min", "ph_mean", "so_mean")  # List all env variables here

# Loop through environmental variables and rasterize each dataset
for (var in env_variables) {
  
  # Ensure that each env variable column exists
  if (var %in% colnames(bio_oracle_2010)) {
    
    # Select the relevant column from bio_oracle_sites_present
    env_data <- bio_oracle_2010 %>% dplyr::select(longitude, latitude, all_of(var))
    
    # Convert the data to matrix of coordinates and values
    coordinates <- cbind(env_data$longitude, env_data$latitude)
    
    # Rasterize the data onto the study_raster template
    env_raster <- rasterize(coordinates, study_raster, env_data[[var]], fun = mean, na.rm = TRUE)
    common_zlim <- c(10, 26)
    
    # Plot the raster for the current environmental variable
    pdf(paste("output/figures/GEA/enviro/Bio-oracle/Raster_bio-oracle_present_",var,".pdf", sep = ""), width = 3.1, height = 5)
    plot(env_raster, xlim = c(-130, 115), ylim = c(33, 46), col = mycolors, axes = TRUE, box = FALSE, xlab="Longitude", ylab="Latitude") 
    dev.off()
 
    # Save the raster to a file
    raster_filename <- paste0("data/processed/GEA/enviro_data/Bio-oracle/tif_files/bio-oracle_present_", var, "_raster.tif")
    writeRaster(env_raster, filename = raster_filename, format = "GTiff", overwrite = TRUE)
  } else {
    warning(paste("Environmental variable", var, "not found in bio_oracle_2010"))
  }
}

# ================================================================================== #

# Save full path of tif files as envtif
env_tif <- list.files("data/processed/GEA/enviro_data/Bio-oracle/tif_files/", pattern ="bio-oracle_present_", full.names = TRUE)

# Stack tif files
env_stack <- stack(env_tif)

# ================================================================================== #

# Convert the raster stack to a dataframe
env_values <- as.data.frame(raster::extract(env_stack, 1:ncell(env_stack), df = TRUE))
env_values$cell <- 1:nrow(env_values)

# ================================================================================== #

# Save the raster stack
writeRaster(env_stack, "data/processed/GEA/enviro_data/Bio-oracle/tif_files/bio-oracle_env_present_stack_raster.tif", format = "GTiff")

# ================================================================================== #
# ================================================================================== #

# Create rasters for future data

# Loop through environmental variables and rasterize each dataset
for (var in env_variables) {
  
  # Ensure that each env variable column exists
  if (var %in% colnames(bio_oracle_ssp585_2090)) {
    
    # Select the relevant column from bio_oracle_sites_present
    env_data <- bio_oracle_ssp585_2090 %>% dplyr::select(longitude, latitude, all_of(var))
    
    # Convert the data to matrix of coordinates and values
    coordinates <- cbind(env_data$longitude, env_data$latitude)
    
    # Rasterize the data onto the study_raster template
    env_raster <- rasterize(coordinates, study_raster, env_data[[var]], fun = mean, na.rm = TRUE)
    common_zlim <- c(10, 26)
    
    # Plot the raster for the current environmental variable
    pdf(paste("output/figures/GEA/enviro/Bio-oracle/Raster_bio-oracle_ssp585_",var,".pdf", sep = ""), width = 3.1, height = 5)
    plot(env_raster, xlim = c(-130, 115), ylim = c(33, 46), col = mycolors, axes = TRUE, box = FALSE, xlab="Longitude", ylab="Latitude") 
    dev.off()
 
    # Save the raster to a file
    raster_filename <- paste0("data/processed/GEA/enviro_data/Bio-oracle/tif_files/bio-oracle_ssp585", var, "_raster.tif")
    writeRaster(env_raster, filename = raster_filename, format = "GTiff", overwrite = TRUE)
  } else {
    warning(paste("Environmental variable", var, "not found in bio_oracle_ssp585_2090"))
  }
}

# ================================================================================== #

# Save full path of tif files as envtif
env_ssp585_tif <- list.files("data/processed/GEA/enviro_data/Bio-oracle/tif_files/", pattern ="bio-oracle_ssp585", full.names = TRUE)

# Stack tif files
env_ssp585_stack <- stack(env_ssp585_tif)

# ================================================================================== #

# Convert the raster stack to a dataframe
env_ssp585_values <- as.data.frame(raster::extract(env_ssp585_stack, 1:ncell(env_ssp585_stack), df = TRUE))
env_ssp585_values$cell <- 1:nrow(env_ssp585_values)

# ================================================================================== #

# Save the raster stack
writeRaster(env_ssp585_stack, "data/processed/GEA/enviro_data/Bio-oracle/tif_files/bio-oracle_env_ssp585_stack_raster.tif", format = "GTiff")

# ================================================================================== #

