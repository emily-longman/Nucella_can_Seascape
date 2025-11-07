# Summarize mussel shell thickness data

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
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
# ================================================================================== #

# Load data

# Read in mussel shell data
Mcali <- read.csv("data/raw/M.cali.thk/Mcali_shell_thk.csv", header=T)

# ================================================================================== #

# Create bin column
Mcali <- Mcali %>% mutate(
        value_bin = case_when(
          Length.L..mm. <= 20 ~ "10",
          Length.L..mm. > 20 & Length.L..mm. <= 30 ~ "20",
          Length.L..mm. > 30 & Length.L..mm. <= 40 ~ "30",
          Length.L..mm. > 40 & Length.L..mm. <= 50 ~ "40",
          Length.L..mm. > 50 & Length.L..mm. <= 60 ~ "50",
          Length.L..mm. > 60 & Length.L..mm. <= 70 ~ "60",
          Length.L..mm. > 70 & Length.L..mm. <= 80 ~ "70",
          Length.L..mm. > 80 & Length.L..mm. <= 90 ~ "80",
          Length.L..mm. > 90 & Length.L..mm. <= 100 ~ "90",
          Length.L..mm. > 100 & Length.L..mm. <= 110 ~ "100",
          Length.L..mm. > 110 & Length.L..mm. <= 120 ~ "110",
          Length.L..mm. > 120 & Length.L..mm. <= 130 ~ "120",
          Length.L..mm. > 130 & Length.L..mm. <= 140 ~ "130",
          Length.L..mm. > 140 & Length.L..mm. <= 150 ~ "140",
          Length.L..mm. > 150 ~ "150",
          TRUE ~ "Other" # Catch-all for values not fitting previous conditions
        )
      )

# ================================================================================== #

# Summarize for total number of mussel collected
Mcali.sum.collected <- Mcali %>% group_by(Site.Code, value_bin) %>% summarize(count=n())

# Write table
write.csv(Mcali.sum.collected, "data/raw/M.cali.thk/Mcali.sum.collected.csv", row.names=FALSE)

# ================================================================================== #

# Filter for only mussels that have done ImageJ on
Mcali.data <- Mcali %>% filter(!is.na(Segment.Area))

# Summarize for total number of mussel collected
Mcali.data.sum <- Mcali.data %>% group_by(Site.Code, value_bin) %>% summarize(count=n())

# Write table
write.csv(Mcali.data.sum, "data/raw/M.cali.thk/Mcali.data.sum.csv", row.names=FALSE)

