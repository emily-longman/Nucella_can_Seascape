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
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]
mycolors_sub5 <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6,7)]

# ================================================================================== #

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# ================================================================================== #

# Load data and merge files - SSP1-1.9

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_119/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_119/per_simcycle/', x))))

# Read all the files and perform ABC
future119_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)
    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_119/per_simcycle/phclineAFs_future_freq.")) %>% 
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

# Join with ecovars
future119_sim_data <- left_join(future119_sim_data, ecovars, by = "sim_eq")

# Save output
save(future119_sim_data, file = "data/processed/SLiM/ph_future/sim_future119_data.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future119_data.Rdata")

#####

# Load data and merge files - SSP1-2.6

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_126/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_126/per_simcycle/', x))))

# Read all the files and perform ABC
future126_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)
    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_126/per_simcycle/phclineAFs_future_freq.")) %>% 
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

# Join with ecovars
future126_sim_data <- left_join(future126_sim_data, ecovars, by = "sim_eq")

# Save output
save(future126_sim_data, file = "data/processed/SLiM/ph_future/sim_future126_data.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future126_data.Rdata")

#####

# Load data and merge files - SSP2-4.5

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_245/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_245/per_simcycle/', x))))

# Read all the files and perform ABC
future245_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)
    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_245/per_simcycle/phclineAFs_future_freq.")) %>% 
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

# Join with ecovars
future245_sim_data <- left_join(future245_sim_data, ecovars, by = "sim_eq")

# Save output
save(future245_sim_data, file = "data/processed/SLiM/ph_future/sim_future245_data.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future245_data.Rdata")

#####
################################
# Load data and merge files - SSP3-7.0

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_370/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_370/per_simcycle/', x))))

# Read all the files and perform ABC
future245_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)
    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_245/per_simcycle/phclineAFs_future_freq.")) %>% 
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

# Join with ecovars
future245_sim_data <- left_join(future245_sim_data, ecovars, by = "sim_eq")

# Save output
save(future245_sim_data, file = "data/processed/SLiM/ph_future/sim_future245_data.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future245_data.Rdata")



#####

# Load data and merge files - SSP5-8.5

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_585/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_585/per_simcycle/', x))))

# Read all the files and perform ABC
future585_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)
    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_585/per_simcycle/phclineAFs_future_freq.")) %>% 
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

# Join with ecovars
future585_sim_data <- left_join(future585_sim_data, ecovars, by = "sim_eq")

# Save output
save(future585_sim_data, file = "data/processed/SLiM/ph_future/sim_future585_data.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future585_data.Rdata")

# ================================================================================== #

# Graph AFs through time

# Make Site factor
future_sim_data$Site <- factor(future_sim_data$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph all sim AFs through time
pdf("output/figures/SLiM/ph_future/AF_time_all.pdf", width = 11.5, height = 7)
ggplot(future_sim_data, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Subset
future_sim_data_sub <- future_sim_data %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")

# Graph subset
pdf("output/figures/SLiM/ph_future/AF_time_all_sub_5pop.pdf", width = 11.5, height = 7)
ggplot(future_sim_data_sub, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Average across the iterations
future_avg <- future_sim_data %>% group_by(year, sim_eq) %>% reframe(AF_fut_avg = mean(AF_fut))

# Join
future_avg <- left_join(ecovars, future_avg, by = "sim_eq")

# Make Site factor
future_avg$Site <- factor(future_avg$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Graph mean (of reps) AFs through time

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time.pdf", width = 10, height = 8.5)
ggplot(future_avg, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    scale_color_manual(values = mycolors) + scale_y_continuous(limits=c(0,1.0), breaks = c(0, 0.5, 1)) +
    labs(x = "Year", y = "Allele Frequency") + theme_linedraw(base_size = 32) + theme(legend.position = "none")
dev.off()
# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_wider.pdf", width = 12, height = 7)
ggplot(future_avg, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    scale_color_manual(values = mycolors) + scale_y_continuous(limits=c(0,1.0), breaks = c(0, 0.5, 1)) + 
    labs(x = "Year", y = "Allele Frequency") + theme_linedraw(base_size = 32) + theme(legend.position = "none")
dev.off()

# Subset
future_avg_sub <- future_avg %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_sub.pdf", width = 9, height = 8)
ggplot(future_avg_sub, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    scale_color_manual(values = mycolors_sub5) + scale_y_continuous(limits=c(0,1.0), breaks = c(0, 0.5, 1)) +
    labs(x = "Year", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Delta AF

# Final AF
final_AF <- future_avg$AF_fut_avg[which(future_avg$year == 100)]
# Add metadata
ph_AF_change <- data.table(ecovars, final_AF)

# Calc delta AF
ph_AF_change$delta_AF <- future_avg$AF_fut_avg[which(future_avg$year == 100)] - future_avg$AF_fut_avg[which(future_avg$year == 21)]

# Make Site factor
ph_AF_change$Site <- factor(ph_AF_change$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph delta AF as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph Delta AF
pdf("output/figures/SLiM/ph_future/Delta_AF_map.pdf", width = 9, height = 9)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = ph_AF_change, aes(x = Longitude, y = Latitude, fill = delta_AF), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Purples"), name="delta AF", breaks = c(0.0, 0.3, 0.6, 0.9)) +
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
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 20), legend.position = c(0.75, 0.825), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
