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
bio_oracle <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_larger_area.csv", header=T)
bio_oracle_future <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_future_larger_area.csv", header=T)

# Extract the most recent data (i.e., 2010-01-01T00:00:00Z, which represents 2010-~2020)
bio_oracle_2010 <- bio_oracle %>% 
  filter(time == "2010-01-01T00:00:00Z")

bio_oracle_future <- bio_oracle_future %>% 
  filter(time == "2090-01-01T00:00:00Z")

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(32, 46.5)
longitude_range <- c(-126, -116)

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
pdf("output/figures/enviro/Bio-oracle/Raster_template.pdf", width = 5, height = 5)
plot(study_raster, main = "West Coast Raster Template")
dev.off()

# ================================================================================== #

# Test graphing just one environmental variable (e.g., mean temperature)

# Set coordinates of bio-oracle data
coordinates <- cbind(bio_oracle_2010$longitude, bio_oracle_2010$latitude)

# Extract one variable - temperature 
bio_oracle_2010_temp_mean <- bio_oracle_2010 %>% 
  dplyr::select(longitude, latitude, thetao_mean)
# Extract one variable - ph mean 
bio_oracle_2010_ph_mean <- bio_oracle_2010 %>% 
  dplyr::select(longitude, latitude, ph_mean)

# Extract one future variable - ph mean 
bio_oracle_future_ph_mean <- bio_oracle_future %>% 
  dplyr::select(longitude, latitude, ph_mean)

# Rasterize temperature data
temp_raster <- rasterize(coordinates, study_raster, bio_oracle_2010_temp_mean$thetao_mean, fun = mean, na.rm = TRUE)
# Rasterize ph data
ph_raster <- rasterize(coordinates, study_raster, bio_oracle_2010_ph_mean$ph_mean, fun = mean, na.rm = TRUE)

# Rasterize ph data
ph_future_raster <- rasterize(coordinates, study_raster, bio_oracle_future_ph_mean$ph_mean, fun = mean, na.rm = TRUE)

# Look at structure
str(temp_raster)

# Graph temperature raster
pdf("output/figures/enviro/Bio-oracle/Test_Raster_bio-oracle_temp_mean.pdf", width = 3.25, height = 5)
plot(temp_raster, col = mycolors, axes = TRUE, box = FALSE,
xlim = c(-130, 114), ylim = c(33, 46), 
xlab="Longitude", ylab="Latitude", main = "Rasterized thetao_mean")
dev.off()
# Graph ph raster
pdf("output/figures/enviro/Bio-oracle/Test_Raster_bio-oracle_ph_mean.pdf", width = 3.25, height = 5)
plot(ph_raster, col = mycolors, axes = TRUE, box = FALSE,
xlim = c(-125, 114), ylim = c(32, 46), 
xlab="Longitude", ylab="Latitude", main = "Rasterized ph_mean")
dev.off()

# Graph future ph raster
pdf("output/figures/enviro/Bio-oracle/Test_Raster_bio-oracle_ph_mean_future.pdf", width = 3.25, height = 5)
plot(ph_future_raster, col = mycolors, axes = TRUE, box = FALSE,
xlim = c(-125, 114), ylim = c(32, 46), 
xlab="Longitude", ylab="Latitude", main = "Rasterized ph_mean future")
dev.off()

# Write raster tif file
writeRaster(temp_raster, filename = "output/figures/enviro/Bio-oracle/Test_biooracle_thetao_mean_test_raster.tif", format = "GTiff")

# Graph with ggplot
# Change to data frame
raster_df_temp <- as.data.frame(temp_raster, xy = TRUE, na.rm = TRUE)
raster_df_ph <- as.data.frame(ph_raster, xy = TRUE, na.rm = TRUE)
# Future
raster_df_ph_future <- as.data.frame(ph_future_raster, xy = TRUE, na.rm = TRUE)
# Join raster files and calc diff
raster_df_ph.tmp <- raster_df_ph %>% rename(current = layer)
raster_df_ph_future.tmp <- raster_df_ph_future %>% rename(future = layer)
raster_df_ph_diff <- left_join(raster_df_ph.tmp, raster_df_ph_future.tmp)
raster_df_ph_diff <- raster_df_ph_diff %>% mutate(diff = future-current)

# Graph temp
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_temp_mean_ggplot.pdf", width = 3.5, height = 5)
ggplot(raster_df_temp, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_fill_gradientn(colours=brewer.pal(5, "RdBu")) +
  #scale_fill_viridis_c() +  # Color scale for the raster values
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  ggtitle("Rasterized ph_mean")  +
  theme_void() +
  theme(legend.title = element_blank(), plot.title = element_text(hjust=0.5))
dev.off()
# Graph pH (note: reversed colors since low pH is the stressor)
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_ph_mean_ggplot_alt.pdf",  width = 6, height = 8) 
ggplot(raster_df_ph, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="mean pH", breaks = c(7.92, 7.96, 8.0, 8.04)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_bw(base_size = 32) + #xlim(c(-126, -117)) + ylim(c(32, 47)) +
  labs(x = "Longitude", y = "Latitude") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.75, 0.53), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
  #theme(plot.title = element_text(hjust=0.5))
dev.off()


# Graph pH future (note: reversed colors since low pH is the stressor)
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_ph_mean_FUTURE_ggplot_alt.pdf",  width = 6, height = 8) 
ggplot(raster_df_ph_future, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="ph future", breaks = c(7.62, 7.66, 7.70, 7.74)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_bw(base_size = 27) + #xlim(c(-126, -117)) + ylim(c(32, 47)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = c(0.75, 0.53), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
  #theme(plot.title = element_text(hjust=0.5))
dev.off()


# Graph pH future (note: reversed colors since low pH is the stressor)
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_ph_mean_TEMPORAL_DIFF_ggplot_alt.pdf",  width = 6, height = 8) 
ggplot(raster_df_ph_diff, aes(x = x, y = y, fill = diff)) +
  geom_raster(aes(fill=diff)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="ph diff") +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_bw(base_size = 27) + #xlim(c(-126, -117)) + ylim(c(32, 47)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = c(0.75, 0.53), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
  #theme(plot.title = element_text(hjust=0.5))
dev.off()
