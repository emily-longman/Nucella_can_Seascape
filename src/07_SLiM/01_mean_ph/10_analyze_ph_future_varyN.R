# Analyze pH future vary N

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

# Load data and merge files for each N

# N values
Ns <- c(2500, 3750, 5000, 10000)

# Loop through each N and extract files and merge
future_sim_data_varyN <- foreach(N=Ns, .combine="rbind", .errorhandling = "remove")%do%{  

  # Create list of file names for data
  file_names = as.list(dir(path = paste0('data/processed/SLiM/ph_future/ph_vary_N/N_', N ,'/'), pattern = "phclineAFs_future_freq.*"))
  file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste0('data/processed/SLiM/ph_future/ph_vary_N/N_', N ,'/'), x))))

  # Read all the files and perform ABC
  future_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = paste0('data/processed/SLiM/ph_future/ph_vary_N/N_', N ,'/phclineAFs_future_freq.'))) %>% 
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

}

# Save output
save(future_sim_data_varyN, file = "data/processed/SLiM/ph_future/ph_vary_N/sim_future_data_varyN.Rdata")
#load("data/processed/SLiM/ph_future/ph_vary_N/sim_future_data_varyN.Rdata")

# ================================================================================== #

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Join
future_sim_data_varyN %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_varyN$Site <- factor(future_sim_data_varyN$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Average across the iterations
future_avg <- future_sim_data_varyN %>% group_by(N, year, sim_eq) %>% reframe(AF_fut_avg = mean(AF_fut))

# Join
future_avg <- left_join(ecovars, future_avg, by = "sim_eq")
# Make Site factor
future_avg$Site <- factor(future_avg$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Order N
future_avg$N <- factor(future_avg$N, levels = c(2500, 3750, 5000, 10000))

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_varyN.pdf", width = 10, height = 16)
ggplot(future_avg, aes(x = year, y = AF_fut_avg, color = Site)) + geom_line(linewidth = 3) + 
    facet_wrap(~N, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

pdf("output/figures/SLiM/ph_future/AF_time_varyN_all.pdf", width = 10, height = 16)
ggplot(future_sim_data_varyN, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    facet_wrap(~N, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Delta AF and level maladapted

# Current AF
current_AF <- future_sim_data_varyN %>% filter(year == 21) %>% rename(AF_2021 = AF_fut)
# Final AF
final_AF <- future_sim_data_varyN %>% filter(year == 100) %>% rename(AF_2100 = AF_fut)

# Join
ph_AF_change <- left_join(current_AF[,-8], final_AF[,-8])

# Calc delta AF
ph_AF_change %<>% mutate(delta_AF = AF_2100-AF_2021)

# Cal level maladapted (assumes AF=1 maximizes fitness at 2100)
ph_AF_change %<>% mutate(mal = 1-AF_2100)

# Make Site factor
ph_AF_change$Site <- factor(ph_AF_change$Site, levels=rev(c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")))
# Make N factor 
s$N <- factor(ph_AF_change$N, levels = c(2500, 3750, 5000, 10000))

# Graph delta AF
pdf("output/figures/SLiM/ph_future/Delta_AF_varyN.pdf", width = 16, height = 10)
ggplot(ph_AF_change, aes(y = Site, x = delta_AF, color = Site)) + geom_boxplot() +
    facet_wrap(~N, ncol = 4) +
    scale_color_manual(values = rev(mycolors)) + 
    labs(x = "Delta AF", y = "Site") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()