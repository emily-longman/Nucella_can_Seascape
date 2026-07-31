# Analyze pH future vary m

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

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/ph_future")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# ================================================================================== #

# Load data and merge files for each m

# m values
ms <- c(0.01, 0.001, 0.0001)

# Loop through each N and extract files and merge
future_sim_data_vary_m <- foreach(m=ms, .combine="rbind", .errorhandling = "remove")%do%{  

  # Create list of file names for data
  file_names = as.list(dir(path = paste0('data/processed/SLiM/ph_future/ph_vary_m/m_', m ,'/'), pattern = "phclineAFs_future_freq.*"))
  file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste0('data/processed/SLiM/ph_future/ph_vary_m/m_', m ,'/'), x))))

  # Read all the files and perform ABC
  future_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = paste0('data/processed/SLiM/ph_future/ph_vary_m/m_', m ,'/phclineAFs_future_freq.'))) %>% 
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
save(future_sim_data_vary_m, file = "data/processed/SLiM/ph_future/ph_vary_m/sim_future_data_vary_m.Rdata")
#load("data/processed/SLiM/ph_future/ph_vary_m/sim_future_data_vary_m.Rdata")

# ================================================================================== #

# Join
future_sim_data_vary_m %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m$Site <- factor(future_sim_data_vary_m$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Average across the iterations
future_avg <- future_sim_data_vary_m %>% group_by(m, year, sim_eq) %>% reframe(AF_fut_avg = mean(AF_fut))

# Join
future_avg <- left_join(ecovars, future_avg, by = "sim_eq")
# Make Site factor
future_avg$Site <- factor(future_avg$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Graph AFs through time

# Order m
future_avg$m <- factor(future_avg$m, levels = c(0.01, 0.001, 0.0001))

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_vary_m.pdf", width = 10, height = 16)
ggplot(future_avg, aes(x = year, y = AF_fut_avg, color = Site)) + geom_line(linewidth = 3) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Order m
future_sim_data_vary_m$m <- factor(future_sim_data_vary_m$m, levels = c(0.01, 0.001, 0.0001))

# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_all.pdf", width = 10, height = 16)
ggplot(future_sim_data_vary_m, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Delta AF and level maladapted

# Current AF
current_AF <- future_sim_data_vary_m %>% filter(year == 21) %>% rename(AF_2021 = AF_fut)
# Final AF
final_AF <- future_sim_data_vary_m %>% filter(year == 100) %>% rename(AF_2100 = AF_fut)

# Join
ph_AF_change <- left_join(current_AF[,-8], final_AF[,-8])

# Calc delta AF
ph_AF_change %<>% mutate(delta_AF = AF_2100-AF_2021)

# Cal level maladapted (assumes AF=1 maximizes fitness at 2100)
ph_AF_change %<>% mutate(mal = 1-AF_2100)

# Make Site factor
ph_AF_change$Site <- factor(ph_AF_change$Site, levels=rev(c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")))
# Make m factor 
ph_AF_change$m <- factor(ph_AF_change$m, levels = c(0.01, 0.001, 0.0001))

# Graph delta AF
pdf("output/figures/SLiM/ph_future/Delta_AF_vary_m.pdf", width = 16, height = 10)
ggplot(ph_AF_change, aes(y = Site, x = delta_AF, color = Site)) + geom_boxplot() +
    facet_wrap(~m, ncol = 4) +
    scale_color_manual(values = rev(mycolors)) + 
    labs(x = "Delta AF", y = "Site") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# When does the allele even get to the northern most populations

min_func <- function(x) {
if (all(is.na(x))) (100) else min(x, na.rm = TRUE)
}

min_func <- function(x) {
if (all(is.na(x))) (NA) else min(x, na.rm = TRUE)
}

# Loop through the iterations and calc when allele gets to SLR and FC
future_sim_north <-
    foreach(mig=levels(future_sim_data_vary_m$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when SLR and FC first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        #SH = min_func(tmp2$year[which(tmp2$Site == "SH" & tmp2$AF_fut > 0)])+2000,
        SLR = min_func(tmp2$year[which(tmp2$Site == "SLR" & tmp2$AF_fut > 0)])+2000,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000)
    }
}

# Reformat
future_sim_north_melt <- future_sim_north %>%
  reshape2::melt(id = c("repId","m"),
                 variable.name = "Site",
                 value.name = "year")

# Make m factor
future_sim_north_melt$m <- factor(future_sim_north_melt$m, levels = c(0.01, 0.001, 0.0001))

# Summarize
future_sim_north_melt_sum <- future_sim_north_melt %>% group_by(Site, m) %>% summarize(mean = mean(year, na.rm = TRUE), sd = sd(year, na.rm = TRUE))

# Graph
pdf("output/figures/SLiM/ph_future/AF_OR_pops_vary_m.pdf", width = 6.5, height = 8.5)
ggplot(future_sim_north_melt, aes(y = year, x = m, fill = Site, color = Site)) + geom_violin() +
    scale_fill_manual(values = rev(mycolors)[18:19]) + scale_color_manual(values = rev(mycolors)[18:19]) + 
    ylim(2022, 2100) +
    labs(x = "Migration", y = "Year") + theme_linedraw(base_size = 32) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.title = element_blank(), legend.position = c(0.2, 0.905), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# How many sims didn't get the allele by 2100 in SLR and FC?
OR_sum_NA <- future_sim_north_melt %>%
  group_by(Site, m) %>%
  summarise(na_count = sum(is.na(year)), na_perc = na_count/length(year)*100)
