# Download Bio-Oracle data (https://www.bio-oracle.org/index.php )
# Use R package biooracler to access the Bio-Oracle dataset via ERDDAP server (https://github.com/bio-oracle/biooracler)
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
#raw_data_path_from_root <- find_root_file("data", "raw", "Bio-oracle", criterion = has_file("README.md"))
# Set working directory as path from root
#setwd(raw_data_path_from_root)

# ================================================================================== #

# Generate output directories
out_dir_present <- paste("data/raw/Bio-oracle/present")
if (!dir.exists(out_dir_present)) {dir.create(out_dir_present)}

# ================================================================================== #

# Load packages
library(tidyverse)
#install.packages('devtools')
library(devtools) 
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)

# ================================================================================== #

# Example online (https://github.com/bio-oracle/biooracler)

# Extract present day data for range of sites

# Define time, lat, and long (set lat and long to fully encompass all sites)
time = c('2000-01-01T00:00:00Z', '2010-01-01T00:00:00Z')
latitude = c(22, 46)
longitude = c(-108, -130)

# Set constraints
constraints = list(time, latitude, longitude)
names(constraints) = c("time", "latitude", "longitude")

# ================================================================================== #

# Present day data

# Set dataset IDs
pH_Present<-"ph_baseline_2000_2018_depthsurf"
SST_Present <-"thetao_baseline_2000_2019_depthsurf"

# Specify datasets and summary statistics to download 
datasets <- list(
  list(dataset_id = SST_Present, variables = c("thetao_max","thetao_min", "thetao_range", "thetao_mean"),
       constraints = constraints, fmt = "csv", directory = out_dir_present),
  list(dataset_id = pH_Present, variables = c("ph_min", "ph_mean"),
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
    new_name <- file.path(out_dir_present, paste0(dataset_id, "_taller_area.csv"))
    file.rename(new_file, new_name)
    message("Rename ", new_file, " to ", new_name)
  } else {
    message("No new file found for ", dataset_id, " or multiple new files detected.")
  }
}

