# Summarize MARINe data for N. canaliculata collection sites
# Note: prior to running the R script, need to load R module on the VACC (module load R/4.4.1)

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)


# ================================================================================== #

# Read in MARINe Biodiversity data
point_contact <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_point_contact_summary.csv", header=T)
quadrat <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_quadrat_summary.csv", header=T)
swath <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_swath_summary.csv", header=T)
