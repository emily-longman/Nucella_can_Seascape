# Merge SLiM output

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
#install.packages(c('data.table', 'tidyverse', 'foreach'))
library(data.table)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Load data
# Create list of file names for  data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_results/', pattern = "phclineAFs.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_results/', x))))

# Read all the files and add columns associated with the iteration and parameters
ph.out <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o <- read.table(i)
    # Add column with file name identifier
    o <- o %>% mutate(file_name = i)
    # Remove parts of name
    o <- o %>% mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_results/phclineAFs."))
    o <- o %>% mutate(file_name = str_remove(file_name, pattern = ".txt"))
    # Separate columns based on parameters
    o <- separate_wider_delim(o, cols = file_name, delim = "_", names = c("repId", "m", "thresh", "k_1", "k_2", "N", "state", "sim.cycle"))
}

# Save output
save(ph.out, file = "data/processed/SLiM/ph_results.Rdata")

# ================================================================================== #
