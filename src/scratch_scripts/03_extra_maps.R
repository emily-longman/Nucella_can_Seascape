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
library(sf)
library(ggnewscale)
library(colorspace)
library(rnaturalearth)
library(tidyverse)
#install.packages('devtools')
library(devtools) 
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)


# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub5 <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6,7)]

# ================================================================================== #

out_dir_present <- paste("data/raw/Bio-oracle/present")
if (!dir.exists(out_dir_present)) {dir.create(out_dir_present)}

# ================================================================================== #

# Extract present day data for range of sites

# Define time, lat, and long (set lat and long to fully encompass all sites)
time = c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
latitude = c(26, 46)
longitude = c(-83, -65)

# Set constraints
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

# ================================================================================== #

# Set dataset IDs
so_Present<-"so_baseline_2000_2019_depthsurf"
SST_Present <-"thetao_baseline_2000_2019_depthsurf"

# Specify datasets and summary statistics to download 
datasets <- list(
   list(dataset_id = SST_Present, variables = c("thetao_max","thetao_min", "thetao_range", "thetao_mean"),
       constraints = constraints, fmt = "csv", directory = out_dir_present),
  list(dataset_id = so_Present, variables = c("so_mean"),
       constraints = constraints, fmt = "csv", directory = out_dir_present))


# Download datasets using for loop
for (dataset in datasets) {
  
  dataset_id <- dataset$dataset_id
  variables <- dataset$variables
  constraints <- dataset$constraints
  
  # List files before downloading
  files_before <- list.files(out_dir_present, full.names = TRUE)
  # Download dataset
  download_layers(dataset_id, variables = variables, constraints = constraints, fmt = "csv", directory = out_dir_present)
  # List files after downloading
  files_after <- list.files(out_dir_present, full.names = TRUE)
  # Identify the newly downloaded file
  new_file <- setdiff(files_after, files_before)
  
  # Rename the file to match the dataset ID
  if (length(new_file) == 1) {
    new_name <- file.path(out_dir_present, paste0(dataset_id, "_east_coast.csv"))
    file.rename(new_file, new_name)
    message("Rename ", new_file, " to ", new_name)
  } else {
    message("No new file found for ", dataset_id, " or multiple new files detected.")
  }
}

# ================================================================================== #

# Read in Bio-oracle present data
salinity <- read.csv("data/raw/Bio-oracle/present/so_baseline_2000_2019_depthsurf_east_coast.csv", header=T)
SST <- read.csv("data/raw/Bio-oracle/present/thetao_baseline_2000_2019_depthsurf_east_coast.csv", header=T)

# Remove first row (i.e. units)
salinity <- salinity[-1,]
SST <- SST[-1,]

# Change latitude and longitude to numeric
salinity <- salinity %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude), so_mean = as.numeric(so_mean))
SST <- SST %>%
mutate(latitude = as.numeric(latitude), longitude = as.numeric(longitude), thetao_mean = as.numeric(thetao_mean))

# Extract the most recent data (i.e., 2010-01-01T00:00:00Z, which represents 2010-~2020)
salinity_2010 <- salinity %>% 
  filter(time == "2010-01-01T00:00:00Z")
SST_2010 <- SST %>% 
  filter(time == "2010-01-01T00:00:00Z")

# ================================================================================== #

# Specify parameters

# Set geographic constraints
latitude_range <- c(26, 46)
#longitude_range <- c(-83, -65)
longitude_range <- c(-83, -67)

# Set study extent
study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])

# Raster resolution (focal cells of bio-oracle data are at 0.05 degree resolution)
raster_resolution <- 0.05

# Create raster layer object - specify study extent, resolution and coordinate reference system
study_raster <- raster(study_extent, res=raster_resolution, crs="+proj=longlat +datum=WGS84")

# Create empty raster template

# Fill dataset with NA
values(study_raster) <- NA

# Graph empty template
pdf("output/figures/enviro/Bio-oracle/Raster_template_east_coast.pdf", width = 5, height = 5)
plot(study_raster, main = "East Coast Raster Template")
dev.off()

# ================================================================================== #

# Set coordinates of bio-oracle data
coordinates <- cbind(salinity_2010$longitude, salinity_2010$latitude)
coordinates <- cbind(SST_2010$longitude, SST_2010$latitude)

# Rasterize data
so_raster <- rasterize(coordinates, study_raster, salinity_2010$so_mean, fun = mean, na.rm = TRUE)
SST_raster <- rasterize(coordinates, study_raster, SST_2010$thetao_mean, fun = mean, na.rm = TRUE)

# Rasterize data
raster_df_so <- as.data.frame(so_raster, xy = TRUE, na.rm = TRUE)
raster_df_SST <- as.data.frame(SST_raster, xy = TRUE, na.rm = TRUE)

# ================================================================================== #

# Atlantic estuary locations
sites <- data.frame(
  site = c(
    "New Meadows River, West Bath",
    "Wells National Estuarine Research Reserve",
    "Plum Island",
    "Narragansett Bay",
    "Hudson River",
    "Delaware Bay",
    "Chesapeake Bay",
    "VA Chesapeake Bay",
    "Currituck Sound",
    "Rachel Carson Reserve",
    "Masonboro Island Reserve",
    "North Inlet-Winyah Bay Reserve",
    "ACE Basin",
    "Jekyll Island, GA",
    "Indian River Lagoon"
  ),
  lon = c(
    -69.892,
    -70.56988,
    -70.795,
    -71.4167,
    -73.9,
    -75.3,
    -76.22469,
    -76.049,
    -75.789,
    -76.639,
    -77.853,
    -79.17079,
    -80.510,
    -81.47661,
    -80.7
  ),
  lat = c(
    43.87969,
    43.310,
    42.71716,
    41.655,
    40.8,
    39.0,
    38.3825,
    36.91156,
    36.233,
    34.706,
    34.13698,
    33.33648,
    32.520,
    31.0968,
    28.2
  ),
  color = c("E","E","E","C","C","C","C","C","C","C","C","C","C","C","C"), 
  shape = c("C","T","C","C","T","C","C","T","C","T","T","C","C","T","C")
)


# Download map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Atlantic coast extent
xmin <- -83
xmax <- -65
ymin <- 26
ymax <- 46

cape_cod_lat <- 41.8

pdf("output/figures/Blue_crabs.pdf", width = 5, height = 5)
ggplot() +
  geom_sf(data = world, fill = "grey90", color = "grey50", linewidth = 0.2) + 
  geom_hline(yintercept = cape_cod_lat, linetype = "dashed", color = "red", linewidth = 1) +
  geom_point(data = sites, aes(x = lon, y = lat), size = 3, color = "navy") +
  #geom_text(data = sites, aes(x = lon, y = lat, label = site),nudge_x = 0.5,size = 4) + 
  annotate("text",x = -81.5, y = cape_cod_lat + 0.4, label = "Historical northern limit (Cape Cod)", color = "red", hjust = 0, size = 4) +
  coord_sf(xlim = c(xmin, xmax), ylim = c(ymin, ymax), expand = FALSE) + labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid.major = element_line(color = "grey85"), axis.title = element_text(size = 12), axis.text = element_text(size = 10))
dev.off()



# Graph with sites
pdf("output/figures/blue_crab.pdf",  width = 8, height = 11) 
ggplot(raster_df_SST, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_fill_gradientn(colours=rev(brewer.pal(9, "Blues"))[1:7], name="SST") +
  new_scale_fill() + 
  geom_point(data = sites, aes(x = lon, y = lat, fill = color, shape = shape), inherit.aes = FALSE, size = 8) +
  scale_fill_manual(values = c("E" = "#ff8800", "C" = "black"), guide = "none") +
  scale_shape_manual(values = c("C" = 21, "T" = 24), guide = "none") +
  geom_hline(yintercept = cape_cod_lat, linetype = "dashed", color = "#787777", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-80, -75, -70)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 16), legend.position = c(0.87, 0.5), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid")) + guides(fill = "none")
dev.off()

pdf("output/figures/blue_crab_nolegend.pdf",  width = 8, height = 11) 
ggplot(raster_df_SST, aes(x = x, y = y, fill = layer)) +
  geom_raster(aes(fill=layer)) +
  scale_fill_gradientn(colours=rev(brewer.pal(9, "Blues"))[1:7],  guide = "none") +
  new_scale_fill() + 
  geom_point(data = sites, aes(x = lon, y = lat, fill = color, shape = shape), inherit.aes = FALSE, size = 8) +
  scale_fill_manual(values = c("E" = "#ff8800", "C" = "black"), guide = "none") +
  scale_shape_manual(values = c("C" = 21, "T" = 24), guide = "none") +
  geom_hline(yintercept = cape_cod_lat, linetype = "dashed", color = "#787777", linewidth = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = c(-80, -75, -70)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(ratio = 1) +  # Fix aspect ratio so the plot is not distorted
  theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 16), legend.position = c(0.87, 0.5), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid")) + guides(fill = "none")
dev.off()

# ================================================================================== #
