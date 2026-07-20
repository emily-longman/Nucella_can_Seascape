# Merge SLiM output

# Clear memory
rm(list=ls())

# Stop exponential
options(scipen = 999)

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'poolfstat', 'foreach'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(poolfstat)
library(foreach)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/ph_future")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data and merge files

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results/per_simcycle/', x))))

# Read all the files and perform ABC
future_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)

    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results/per_simcycle/phclineAFs_future_freq.")) %>% 
            mutate(file_name = str_remove(file_name, pattern = ".txt"))
    # Rename pops
    names(tmp)[1:19] = paste("p", 0:18, sep ="")
    
    # Separate columns based on parameters
    tmp2 <- separate_wider_delim(tmp, cols = file_name, delim = "_", names = c("repId", "m", "thresh", "k", "mag", "N", "state"))

    # Add year
    tmp2$year <- seq(21:100)

    # Reformat
    #sim_Data_melt <- sim_Data %>%
    #        reshape2::melt(id = c("repId","m","thresh","k_1","k_2","N","state","sim.cycle"), 
    #        variable.name = "sim_eq", value.name = "AF_true")

}

# ================================================================================== #

