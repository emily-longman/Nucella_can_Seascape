# Graph data from Sanford 2003

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'maps', 'mapdata'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(maps)
library(mapdata)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/phenotypic_data")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load Phenotypic data

# Read Sanford 2003 data
drilling <- read.csv("data/raw/phenotypic_data/Sanford_2003.csv", header=T)

# ================================================================================== #

# Order
lat <- c("FC", "SLR", "SH", "ARA", "CBL", "STC", "VD", "BMR", "PGP", "SBR", "PB")
drilling.order <- drilling %>% mutate(Site = factor(Site, levels = lat)) %>% arrange(Site)

# ================================================================================== #

# Color palette
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(11))

# Graph drilling data versus latitude
pdf("output/figures/phenotypic/Drilling_lat.pdf", width = 12, height = 8)
ggplot(data = drilling, aes(x = Lat, y = Mussels_drilled_per_whelk)) +
  geom_point(size = 3, fill="black", shape = 21) + 
  geom_vline(xintercept=36.8007, linetype="solid", color="black") +
  ylab("Drilling") + xlab("Latitude") +
  theme_bw(base_size = 28)
dev.off()



# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))


# Graph map
pdf("output/figures/phenotypic/Sanford2003.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = afs.melt[which(afs.melt$SNP_id == "4922150"),], aes(x = Long, y = Lat, fill = AF), shape = 21, size = 5) + 
  scale_fill_gradient(low = "firebrick", high = "black") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 12) + ggtitle("AF 4922150")
dev.off()

