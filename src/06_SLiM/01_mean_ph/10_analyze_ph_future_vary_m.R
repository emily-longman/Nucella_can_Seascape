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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'poolfstat', 'foreach', 'ggplot2', 'RColorBrewer', 'scales'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(poolfstat)
library(foreach)
library(ggplot2)
library(RColorBrewer)
library(scales)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]
mycolors_sub5 <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6,7)]

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
ms <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m/phclineAFs_future_freq.')) %>% 
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


# Fix scientific notation
future_sim_data_vary_m$m[which(future_sim_data_vary_m$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m$m[which(future_sim_data_vary_m$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m$m <- as.numeric(future_sim_data_vary_m$m)

# Make column with scientific notation
future_sim_data_vary_m <- future_sim_data_vary_m %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m, file = "data/processed/SLiM/ph_future/ph_vary_m/sim_future_data_vary_m.Rdata")
#load("data/processed/SLiM/ph_future/ph_vary_m_old/sim_future_data_vary_m.Rdata")

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
future_avg$m <- factor(future_avg$m, levels = ms)

# Graph AF vs time
pdf("output/figures/SLiM/ph_future/AF_time_vary_m.pdf", width = 10, height = 16)
ggplot(future_avg, aes(x = year, y = AF_fut_avg, color = Site)) + geom_line(linewidth = 3) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Order m
future_sim_data_vary_m$m <- factor(future_sim_data_vary_m$m, levels = ms)

# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_all.pdf", width = 10, height = 16)
ggplot(future_sim_data_vary_m, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Graph 5 pop subset

# Subset
future_avg_sub5 <- future_avg %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
future_sim_data_vary_m_sub5 <- future_sim_data_vary_m %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")

# Graph
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_sub_5pop.pdf", width = 10, height = 16)
ggplot(future_avg_sub5, aes(x = (year+2000), y = AF_fut_avg, color = Site)) + geom_line(linewidth = 4) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + 
    theme(panel.spacing.y = unit(0, "lines")) + 
    theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_all_sub_5pop.pdf", width = 10, height = 16)
ggplot(future_sim_data_vary_m_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Subset to only 3 migrations
future_avg_sub_3m <- future_avg_sub5[which(future_avg_sub5$m == 0.001 | future_avg_sub5$m == 0.0001 | future_avg_sub5$m == 0.00001)]
future_avg_sub_3m$m <- factor(future_avg_sub_3m$m)

# Graph with different linetypes
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_sub_linetype_5pop.pdf", width = 17.5, height = 8.5)
ggplot(future_avg_sub_3m, aes(x = (year+2000), y = AF_fut_avg, color = Site, linetype = m)) + 
    geom_line(linewidth = 5) + 
    scale_color_manual(values = mycolors_sub5) + 
    scale_linetype_manual(values = c("dotted", "dashed", "solid")) +
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 32) + 
    theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))
dev.off()

# ================================================================================== #

# When does the allele even get to the northern most populations

min_func <- function(x) {
if (all(is.na(x))) (NA) else min(x, na.rm = TRUE)
}

# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north <-
    foreach(mig=levels(future_sim_data_vary_m$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000)
    }
}

# Reformat
future_sim_north_melt <- future_sim_north %>%
  reshape2::melt(id = c("repId","m"),
                 variable.name = "Site",
                 value.name = "year")

# Make numeric
future_sim_north_melt$m <- as.numeric(future_sim_north_melt$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_melt_NA <- future_sim_north_melt %>%
  group_by(Site, m) %>%
  summarise(na_count = sum(is.na(year)), na_perc = na_count/length(year)*100)

# Make column with scientific notation
future_sim_north_melt_NA <- future_sim_north_melt_NA %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Make scientific function more presentable
fancy_scientific <- function(l) {
  l <- as.character(l)
  l <- gsub("e", " * 0^", l)
  parse(text = l)
}

# For m = 10^-2 in which most of the sim did get the allele when did it arrive?
mean(future_sim_north_melt$year[which(future_sim_north_melt$m == 0.01)])
sd(future_sim_north_melt$year[which(future_sim_north_melt$m == 0.01)])

# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_NAs.pdf", width = 7.5, height = 8.1)
ggplot(future_sim_north_melt_NA, aes(y = na_perc, x = m_scientific, fill = Site, color = Site)) + 
    geom_point(size = 12) +
    scale_fill_manual(values =  mycolors_sub[1]) + scale_color_manual(values =  mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(-5, 105)+
    labs(x = "Migration", y = "% Simulations Allele\nFailed to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()


# Make m a factor
# Make column with scientific notation
future_sim_north_melt <- future_sim_north_melt %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Graph only those that do make it
future_sim_north_melt_subset <- future_sim_north_melt %>% filter(m %in% ms[1:2])

pdf("output/figures/SLiM/ph_future/AF_OR_pops_vary_m_smaller.pdf", width = 3, height = 4)
ggplot(future_sim_north_melt_subset, aes(y = year, x = m_scientific,  color = Site, fill = Site)) + 
    #geom_jitter(col = "black", alpha = 0.6) + 
    geom_violin() +
    scale_fill_manual(values = mycolors_sub[1]) + 
    scale_color_manual(values = mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(2022, 2100) + 
    labs(x = "Migration", y = "") + theme_linedraw(base_size = 24) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()
