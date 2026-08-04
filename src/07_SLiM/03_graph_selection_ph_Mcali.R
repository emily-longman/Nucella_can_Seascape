# Graph selection coefficient for pH and Mcali

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'RColorBrewer'))
library(tidyverse)
library(data.table)
library(magrittr)
library(RColorBrewer)

# ================================================================================== #

# Load data
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# ================================================================================== #

# Selection
sel <- function(x, z, k) {
  1 / (1 + exp((x - z)/k)) - 0.5
}

# Calc selection coefficient for each var
ecovars$pH <- sel(ecovars$ph_mean, 7.986, 0.10)
ecovars$shell_thickness <- sel(ecovars$shell_thk, 1.91, 42)

# Format
sel <- ecovars[,c(2,9,10)] %>%
  reshape2::melt(id = c("Site"),
                 variable.name = "var",
                 value.name = "sel")


# ================================================================================== #

# Graph selection coefficients for the 19 pops
pdf("output/figures/SLiM/selection_ph_Mcali.pdf", width = 5.5, height = 8.295)
ggplot(sel, aes(x = var, y = sel)) + geom_boxplot() + geom_jitter() +
    scale_x_discrete(labels = c("pH" = "Mean pH", "shell_thickness" = "Shell\nThickness")) +
    labs(x= "", y= "Selection coefficient") + theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()