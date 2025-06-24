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
# install.packages(c('data.table', 'tidyverse', 'plyr, 'foreach', 'ggplot2'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)
library(ggplot2)

# ================================================================================== #

# Load data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1/', pattern = "GLM_100perm_Bio-Oracle_*"))
file_names = lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1/', x))

out = lapply(file_names, function(x){
  env = new.env()
  nm = load(x, envir = env)[1]
  objname = gsub(pattern = 'data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1/', replacement = '', x = x, fixed = T)
  objname = gsub(pattern = 'GLM_100perm_Bio-Oracle_|.RData', replacement = '', x = objname)
  assign(objname, env[[nm]], envir = .GlobalEnv)
  0 # succeeded
} )




load("data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1/GLM_100perm_Bio-Oracle_chunk_1_chr_Backbone_17064_start_1568_stop_28650.Rdata")

load("data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1/GLM_100perm_Bio-Oracle_chunk_1_chr_Backbone_17065_start_678_stop_50678.Rdata")

#### Struggling to load each and merge them together

my_files <- list.files(path = "data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_1", pattern = "GLM_100perm_Bio-Oracle_*", full.names = T)

temp <- new.env()
list <- lapply(my_files, load, temp)

all_data <- as.list(temp)
rm(temp) 

glm.model.output.all <- bind_rows(all_data)



str(glm.model.output.all)