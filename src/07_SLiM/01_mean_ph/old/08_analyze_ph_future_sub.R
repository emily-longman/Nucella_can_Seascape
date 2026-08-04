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
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]

# ================================================================================== #

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars_sub <- ecovars %>% filter(Site == "FC" | Site == "CBL" | Site == "BMR" | Site == "STR")
ecovars_sub %<>% mutate(sim_eq = paste("p", 0:3, sep =""))

# ================================================================================== #

# Load data and merge files

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_results_sub/per_simcycle/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_results_sub/per_simcycle/', x))))

# Read all the files and perform ABC
future_sim_data_sub <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which file loading
    message(i)

    # Load file and add file name identifier
    tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_future/ph_results_sub/per_simcycle/phclineAFs_future_freq.")) %>% 
            mutate(file_name = str_remove(file_name, pattern = ".txt"))
    # Rename pops
    names(tmp)[1:4] = paste("p", 0:3, sep ="")
    
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
save(future_sim_data_sub, file = "data/processed/SLiM/ph_future/ph_results/sim_future_data_sub.Rdata")
#load("data/processed/SLiM/ph_future/ph_results/sim_future_data_sub.Rdata")

# Join with ecovars
future_sim_data_sub <- left_join(future_sim_data_sub, ecovars, by = "sim_eq")

# ================================================================================== #

# Graph AFs through time

# Make Site factor
future_sim_data_sub$Site <- factor(future_sim_data_sub$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph all sim AFs through time
pdf("output/figures/SLiM/ph_future/AF_time_all_sub.pdf", width = 12, height = 7)
ggplot(future_sim_data_sub, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    scale_color_manual(values = mycolors_sub) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 32) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Average across the iterations
future_sub_avg <- future_sim_data_sub %>% group_by(year, sim_eq) %>% reframe(AF_fut_avg = mean(AF_fut))

# Join
future_sub_avg <- left_join(ecovars_sub, future_sub_avg, by = "sim_eq")

# Make Site factor
future_sub_avg$Site <- factor(future_sub_avg$Site, levels=c("FC", "CBL", "BMR", "STR"))

# ================================================================================== #

# Graph mean (of reps) AFs through time

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_sub.pdf", width = 9.75, height = 8.5)
ggplot(future_sub_avg, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    scale_color_manual(values = mycolors_sub) + scale_y_continuous(limits=c(0,1.0), breaks = c(0, 0.5, 1)) +
    labs(x = "Year", y = "Allele Frequency") + theme_linedraw(base_size = 32) + theme(legend.position = "none")
dev.off()

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_sub_wider.pdf", width = 12, height = 7)
ggplot(future_sub_avg, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    scale_color_manual(values = mycolors_sub) + scale_y_continuous(limits=c(0,1.0), breaks = c(0, 0.5, 1)) + 
    labs(x = "Year", y = "Allele Frequency") + theme_linedraw(base_size = 32) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Delta AF and level maladapted

# Final AF
final_AF <- future_avg$AF_fut_avg[which(future_avg$year == 100)]
# Add metadata
ph_AF_change <- data.table(ecovars, final_AF)

# Calc delta AF
ph_AF_change$delta_AF <- future_avg$AF_fut_avg[which(future_avg$year == 100)] - future_avg$AF_fut_avg[which(future_avg$year == 21)]

# Cal level maladapted (assumes AF=1 maximizes fitness at 2100)
ph_AF_change %<>% mutate(mal = 1-final_AF)

# Make Site factor
ph_AF_change$Site <- factor(ph_AF_change$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

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


# Graph level of maladapted
pdf("output/figures/SLiM/ph_future/Maladapt_map.pdf", width = 9, height = 9)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = ph_AF_change, aes(x = Longitude, y = Latitude, fill = mal), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greens"), name="Maladapted", breaks = c(0.0, 0.05, 0.10, 0.15)) +
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
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 20), legend.position = c(0.7, 0.825), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()


# ================================================================================== #

# Genetic load

# Selection
s <- function(x, z, k) {
  1 / (1 + exp((x - z)/k)) - (1/2)
}

# Genetic Load Currently

# Allele frequency data for top ph hits
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit and join with eco vars
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595")  %>%
  left_join(dplyr::select(ecovars, Site, Latitude, sim_eq)) %>%
  arrange(-Latitude)

# Loop through pops - QUESTION: unsure if wii should be 1!!!
ph_w_current <- foreach(i=1:19, .combine="rbind", .errorhandling = "remove")%do%{  
  #message(ecovars$Site[i])
  # Extract current pH
  env_i = ecovars$ph_mean[i]
  # Calculate selection
  s_i = s(env_i, 7.996, 0.09)
  # Calc p and q
  p_i = topsnp$AF[i]
  q_i = 1-p_i
  # Calc fitness of pop
  w_i = (p_i)^2*1 + 2*p_i*q_i*(1-0.5*s_i) + (q_i)^2*(1-s_i)

  # Format
  data.frame(
    Site = ph_AF_change$Site[i],
    Latitude = ph_AF_change$Latitude[i],
    Longitude = ph_AF_change$Longitude[i],
    s_i = s_i,
    w = w_i)
}

# Calc L (genetic load)
ph_w_current %<>% mutate(L = (max(ph_w_current$w)-w)/max(ph_w_current$w))

# Graph level of genetic load
pdf("output/figures/SLiM/ph_future/Genetic_load_current_map.pdf", width = 9, height = 9)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = ph_w_current, aes(x = Longitude, y = Latitude, fill = L), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Blues"), name="Genetic Load", breaks = c(0.01, 0.04, 0.07)) +
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
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 20), legend.position = c(0.68, 0.825), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()


######

# Genetic Load in 2100

# Load future ph lms
s <- read.csv("data/processed/SLiM/ph_future/ph_future_lm.csv")

# Loop through pops
ph_w_future <- foreach(i=1:19, .combine="rbind", .errorhandling = "remove")%do%{  
  message(ph_future_lm$Site[i])
  # Extract env at 2100
  env_i = ph_future_lm$slope[i]*80 + ph_future_lm$intercept[i]
  # Calculate selection
  s_i = s(env_i, 7.987, 0.03)
  # Calc p and q
  p_i = ph_AF_change$final_AF[i]
  q_i = 1-p_i
  # Calc fitness of pop
  w_i = (p_i)^2*1 + 2*p_i*q_i*(1-0.5*s_i) + (q_i)^2*(1-s_i)

  # Format
  data.frame(
    Site = ph_AF_change$Site[i],
    Latitude = ph_AF_change$Latitude[i],
    Longitude = ph_AF_change$Longitude[i],
    s_i = s_i,
    w = w_i)
}

# Calc L (genetic load)
ph_w_future %<>% mutate(L = (max(ph_w_future$w)-w)/max(ph_w_future$w))

# Graph level of genetic load
pdf("output/figures/SLiM/ph_future/Genetic_load_future_map.pdf", width = 9, height = 9)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = ph_w_future, aes(x = Longitude, y = Latitude, fill = L), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Blues"), name="Genetic Load", breaks = c(0.01, 0.04, 0.07)) +
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
  theme(legend.title = element_text(size = 22), legend.text = element_text(size = 20), legend.position = c(0.68, 0.825), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Compare gGO with change in AF

# Load gGO
load("data/processed/genomic_offset/Nucella_gGO.Rdata")

# Join
ph_AF_change_go <- left_join(ph_AF_change, go.scaled.output[,c(1,5)], by = "Site")

# Make Site factor
ph_AF_change_go$Site <- factor(ph_AF_change_go$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph delta AF and gGO
pdf("output/figures/SLiM/ph_future/Delta_AF_gGO.pdf", width = 5.5, height = 5)
ggplot(ph_AF_change_go, aes(x = delta_AF, y = GO.scaled, fill = Site)) + geom_point(size = 6, shape = 21, alpha = 0.9) + scale_fill_manual(values = mycolors) +
    scale_x_continuous(limits = c(-0.05, 1), breaks = c(0, 0.5, 1.0)) + scale_y_continuous(limits = c(0.0801, 0.0991), breaks = c(0.08, 0.09, 0.10)) + 
    labs(x = "Delta AF", y = "gGO") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()
