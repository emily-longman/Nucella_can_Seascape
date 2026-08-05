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

  #file_names = as.list(dir(path = paste0('data/processed/SLiM/ph_future/ph_vary_m/m_', m ,'/'), pattern = "phclineAFs_future_freq.*"))
  #file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste0('data/processed/SLiM/ph_future/ph_vary_m/m_', m ,'/'), x))))

  # Read all the files and perform ABC
  #future_sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    
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

# Graph 4 pop subset

# Subset
future_avg_sub <- future_avg %>% filter(Site == "FC" | Site == "CBL" | Site == "BMR" | Site == "STR")
future_sim_data_vary_m_sub <- future_sim_data_vary_m %>% filter(Site == "FC" | Site == "CBL" | Site == "BMR" | Site == "STR")

# Graph
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_sub.pdf", width = 10, height = 16)
ggplot(future_avg_sub, aes(x = year, y = AF_fut_avg, color = Site)) + geom_line(linewidth = 3) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Subset to only 3 migrations
future_avg_sub_test <- future_avg_sub[which(future_avg_sub$m == 0.01 | future_avg_sub$m == 0.001 | future_avg_sub$m == 0.0001)]
future_avg_sub_test$m <- factor(future_avg_sub_test$m)
# Graph with different linetypes
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_sub_linetype.pdf", width = 18, height = 8.5)
ggplot(future_avg_sub_test, aes(x = year, y = AF_fut_avg, color = Site, linetype = m)) + 
    geom_line(linewidth = 5) + 
    scale_color_manual(values = mycolors_sub) + 
    scale_linetype_manual(values = c("dashed", "solid", "dotted")) +
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 32) + 
    theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))
dev.off()

# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_all_sub.pdf", width = 10, height = 16)
ggplot(future_sim_data_vary_m_sub, aes(x = year, y = AF_fut, color = Site, group=interaction(repId, Site))) + geom_line(linewidth = 1, alpha = 0.4) + 
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub) + 
    labs(x = "Years", y = "AF") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()


# ================================================================================== #

# When does the allele even get to the northern most populations

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

# How many sims didn't get the allele by 2100 in SLR and FC?
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

# Graph - for those that got the allele when did it occur
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_NAs.pdf", width = 6.5, height = 5)
ggplot(future_sim_north_melt_NA, aes(y = na_perc, x = m_scientific, fill = Site, color = Site)) + 
    geom_point(size = 8) +
    scale_fill_manual(values =  mycolors_sub[1]) + scale_color_manual(values =  mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(-5, 105)+
    labs(x = "Migration", y = "Percent") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_NAs_taller.pdf", width = 7, height = 8.5)
ggplot(future_sim_north_melt_NA, aes(y = na_perc, x = m_scientific, fill = Site, color = Site)) + 
    geom_point(size = 12) +
    scale_fill_manual(values =  mycolors_sub[1]) + scale_color_manual(values =  mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(-5, 105)+
    labs(x = "Migration", y = "Percent") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()


# Make m a factor
#future_sim_north_melt$m <- factor(future_sim_north_melt$m, levels = rev(ms))
# Make column with scientific notation
future_sim_north_melt <- future_sim_north_melt %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Graph - for those that got the allele when did it occur
pdf("output/figures/SLiM/ph_future/AF_OR_pops_vary_m_shorter.pdf", width = 6.5, height = 5)
ggplot(future_sim_north_melt, aes(y = year, x = m_scientific,  color = Site, fill = Site)) + 
    #geom_jitter(col = "black", alpha = 0.6) + 
    geom_violin() +
    scale_fill_manual(values = mycolors_sub[1]) + 
    scale_color_manual(values = mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(2022, 2100) + 
    labs(x = "Migration", y = "Year") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

pdf("output/figures/SLiM/ph_future/AF_OR_pops_vary_m_taller.pdf", width = 6.5, height = 8.5)
ggplot(future_sim_north_melt, aes(y = year, x = m_scientific,  color = Site, fill = Site)) + 
    #geom_jitter(col = "black", alpha = 0.6) + 
    geom_violin() +
    scale_fill_manual(values = mycolors_sub[1]) + 
    scale_color_manual(values = mycolors_sub[1]) + 
    scale_x_discrete(labels = fancy_scientific) + ylim(2022, 2100) + 
    labs(x = "Migration", y = "Year") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph only those that do make it
future_sim_north_melt_subset <- future_sim_north_melt %>% filter(m %in% ms[1:3])

pdf("output/figures/SLiM/ph_future/AF_OR_pops_vary_m_smaller.pdf", width = 3.2, height = 4)
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








# ================================================================================== #
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
ph_AF_change$m <- factor(ph_AF_change$m, levels = ms)

# Graph delta AF
pdf("output/figures/SLiM/ph_future/Delta_AF_vary_m.pdf", width = 16, height = 10)
ggplot(ph_AF_change, aes(y = Site, x = delta_AF, color = Site)) + geom_boxplot() +
    facet_wrap(~m, ncol = 4) +
    scale_color_manual(values = rev(mycolors)) + 
    labs(x = "Delta AF", y = "Site") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()
