# Graph Bio-Oracle data (https://www.bio-oracle.org/index.php)
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

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
#install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'terra', 'raster', 'ggnewscale', 'colorspace'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(terra)
library(raster)
library(ggnewscale)
library(colorspace)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]
mycolors_sub5 <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6,7)]

# ================================================================================== #

# Figure directory
out_dir_fig <- paste("output/figures/enviro/Bio-oracle")
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

# Metadata
meta <- fread("guide_files/Nucella_ph_shellt.txt")
names(meta)[2] = "Site"

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(32, 46.5)
longitude_range <- c(-126, -118)

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

# Graph empty template
pdf("output/figures/enviro/Bio-oracle/Raster_template.pdf", width = 5, height = 5)
plot(study_raster, main = "West Coast Raster Template")
dev.off()

# ================================================================================== #

# Test graphing just one environmental variable (e.g., mean temperature)

# Set coordinates of bio-oracle data
coordinates <- cbind(bio_oracle_2010$longitude, bio_oracle_2010$latitude)

# Extract one variable - ph mean 
bio_oracle_2010_ph_mean <- bio_oracle_2010 %>% 
  dplyr::select(longitude, latitude, ph_mean)

# Extract one future variable - ph mean 
bio_oracle_future_ph_mean <- bio_oracle_future %>% 
  dplyr::select(longitude, latitude, ph_mean)

# Rasterize ph data
ph_raster <- rasterize(coordinates, study_raster, bio_oracle_2010_ph_mean$ph_mean, fun = mean, na.rm = TRUE)

# Rasterize ph data
ph_future_raster <- rasterize(coordinates, study_raster, bio_oracle_future_ph_mean$ph_mean, fun = mean, na.rm = TRUE)

# Graph with ggplot
# Change to data frame
raster_df_ph <- as.data.frame(ph_raster, xy = TRUE, na.rm = TRUE)
# Future
raster_df_ph_future <- as.data.frame(ph_future_raster, xy = TRUE, na.rm = TRUE)
# Join raster files and calc diff
raster_df_ph.tmp <- raster_df_ph %>% rename(current = layer)
raster_df_ph_future.tmp <- raster_df_ph_future %>% rename(future = layer)
raster_df_ph_diff <- left_join(raster_df_ph.tmp, raster_df_ph_future.tmp)
raster_df_ph_diff <- raster_df_ph_diff %>% mutate(diff = future-current)

# ================================================================================== #

# Graph Delta pH future (note: reversed colors since low pH is the stressor)
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_DELTA_ph_mean_ggplot.pdf",  width = 6, height = 8.4) 
ggplot(raster_df_ph_diff, aes(x = x, y = y, fill = diff)) +
  geom_raster(aes(fill=diff)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="Delta pH") +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.7, 0.82), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
  #theme(plot.title = element_text(hjust=0.5))
dev.off()

# ================================================================================== #

# Graph delta ph with sites

# Subset sites
sites5 <- meta %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")
sites5$Site <- factor(sites5$Site, levels = c("FC", "ARA", "CBL", "BMR", "STR"))

# Graph with sites
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_DELTA_ph_mean_ggplot_sites_sub_alt_larger_5pop.pdf",  width = 7.5, height = 10) 
ggplot(raster_df_ph_diff, aes(x = x, y = y, fill = diff)) +
  geom_raster(aes(fill=diff)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="Delta pH") +
  new_scale_fill() + 
  geom_point(data = sites5, aes(x = Longitude, y = Latitude, fill = Site), shape=21, inherit.aes = FALSE, size = 14) +
  scale_fill_manual(values = mycolors_sub5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 26), legend.text = element_text(size = 20), legend.position = c(0.725, 0.83), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid")) + guides(fill = "none")
dev.off()
pdf("output/figures/enviro/Bio-oracle/Raster_bio-oracle_DELTA_ph_mean_ggplot_sites_sub_alt_smaller_5pop.pdf",  width = 6, height = 8) 
ggplot(raster_df_ph_diff, aes(x = x, y = y, fill = diff)) +
  geom_raster(aes(fill=diff)) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "YlOrRd")), name="Delta pH") +
  new_scale_fill() + 
  geom_point(data = sites5, aes(x = Longitude, y = Latitude, fill = Site), shape=21, inherit.aes = FALSE, size = 11) +
  scale_fill_manual(values = mycolors_sub5) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 26), legend.text = element_text(size = 20), legend.position = c(0.725, 0.82), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid")) + guides(fill = "none")
dev.off()