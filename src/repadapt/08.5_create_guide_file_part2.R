# Create guide file for mpileup.

# ================================================================================== #

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo

# Install and load rprojroot
#install.packages(c('rprojroot'))
library(rprojroot)
# Create relative path from root
rel_path_from_root <- find_root_file("data", "processed", "repadapt", "guide_files" ,criterion = has_file("README.md"))
# Set working directory as path from root
setwd(rel_path_from_root)

# ================================================================================== #

# Load packages
#install.packages(c("data.table", "tidyverse", "groupdata2"))
library(data.table)
library(tidyverse)
library(groupdata2)

# ================================================================================== #

# Read in data
data <- fread("chromosomes.txt", header = F)

# Group the backbone names into 30
group(data, n=30, method = "greedy") -> guide_file_array

# Write the table
write.table(guide_file_array, "chromosomes_guide_file.txt", col.names = F, row.names = F, quote = F)
# Note guide_file_array has dimensions:  19014, 2
