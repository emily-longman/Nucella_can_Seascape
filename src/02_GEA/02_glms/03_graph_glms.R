# Graph glms output

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
install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach', 'ggplot2'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)
library(ggplot2)

# ================================================================================== #

# Load data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_window_analysis/', pattern = "GLM_100perm_Bio-Oracle_chunk_*"))
file_names = lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_window_analysis/', x))
  
# Read all the files 
collated.files =  
foreach(i=file_names, .combine="rbind")%do%{  
message(i) 
chunk = i
o = get(load(i))
o %>% mutate(chunk = file.name) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_"))
} 

