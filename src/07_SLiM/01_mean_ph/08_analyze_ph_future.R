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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'poolfstat', 'foreach', 'ggplot2', 'RColorBrewer'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(poolfstat)
library(foreach)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/ph_future")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

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
    tmp2$year <- c(21:100)

    # Reformat
    tmp2_melt <- tmp2 %>%
            reshape2::melt(id = c("repId","m","thresh","k","mag","N","state","year"), 
            variable.name = "sim_eq", value.name = "AF_fut")
}

# Save output
save(future_sim_data, file = "data/processed/SLiM/ph_future/ph_results/sim_future_data.Rdata")
write.csv(future_sim_data, file = "data/processed/SLiM/ph_future/ph_results/sim_future_data.csv", row.names=F)
#load("data/processed/SLiM/ph_future/ph_results/sim_future_data.Rdata")

# ================================================================================== #

# Average across the iterations
future_avg <- future_sim_data %>% group_by(year, sim_eq) %>% reframe(AF_fut_avg = mean(AF_fut))

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Join
future_avg <- left_join(ecovars, future_avg, by = "sim_eq")

# Make Site factor
future_avg$Site <- factor(future_avg$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Calc diff
delta_AF <- future_avg$AF_fut_avg[which(future_avg$year == 80)] - future_avg$AF_fut_avg[which(future_avg$year == 1)]

# Add metadata
future_diff <- data.table(ecovars, delta_AF)

# Make Site factor
future_diff$Site <- factor(future_diff$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Graph delta AF as map

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph Delta AF
pdf("output/figures/SLiM/ph_future/Delta_AF_map.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = future_diff, aes(x = Longitude, y = Latitude, fill = delta_AF), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"), name="delta AF") +
  coord_fixed(1.3) +
  scale_x_continuous(
    limits = c(-125, -114.1),
    breaks = seq(-125, -114.1, by = 3) # Tick marks every 0.5 units
  ) + ylim(32, 48) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 32) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 14), legend.position = c(0.818, 0.51))
dev.off()

# ================================================================================== #

# TRYING TO ANIMATE OVER TIME BUT NOT CURRENTLY WORKING

# Must have sf so need to load gdal
# module load gdal/3.11.4
# install.packages("sf")
library(sf)
#devtools::install_github('thomasp85/gganimate')
library(gganimate)
library(gifski)

# Show change in AF over time
pdf("output/figures/SLiM/ph_future/AF_time.pdf", width = 8, height = 5)
ggplot(future_avg, aes(x = year, y = AF_fut_avg, group = Site, color = Site)) +
    geom_line(linewidth=2) + scale_color_manual(values = mycolors) +
    labs(x = "Year", y = "AF") + 
    theme_linedraw(base_size = 30) + theme(legend.position="none")
dev.off()

test <- ggplot(future_avg, aes(x = year, y = AF_fut_avg, group = Site, color = Site)) +
    geom_line(linewidth=2) + scale_color_manual(values = mycolors) +
    labs(x = "Year", y = "AF") + 
    theme_linedraw(base_size = 30) + theme(legend.position="none") +  transition_reveal(year)

animate(test, duration = 5, fps = 20, width = 200, height = 200, renderer = file_renderer("output/figures/SLiM/ph_future/AF_test.gif"))

anim_save(animation = test, filename = "output/figures/SLiM/ph_future/AF_test.gif")



# Graph AF changing over time
AF_map <- ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = future_avg, aes(x = Longitude, y = Latitude, fill = AF_fut_avg), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"), name="delta AF") +
  coord_fixed(1.3) +
  scale_x_continuous(
    limits = c(-125, -114.1),
    breaks = seq(-125, -114.1, by = 3) # Tick marks every 0.5 units
  ) + ylim(32, 48) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 32) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 14), legend.position = c(0.818, 0.51)) +
  transition_time(year)

# Save as gif
anim_save(AF_map, "output/figures/processed/genomic_offset/AF_future.gif")



# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Compare gGO with change in AF

# Load gGO
load("data/processed/genomic_offset/Nucella_gGO.Rdata")

# Join
future_diff_go <- left_join(future_diff, go.scaled.output[,c(1,5)], by = "Site")

# Make Site factor
future_diff_go$Site <- factor(future_diff_go$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph delta AF and gGO
pdf("output/figures/SLiM/ph_future/Delta_AF_gGO.pdf", width = 6, height = 5)
ggplot(future_diff_go, aes(x = delta_AF, y = GO.scaled, fill = Site)) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = mycolors) +
    theme_linedraw()
dev.off()





