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

# m values
ms <- c(0.01, 0.001, 0.0001, 0.00001, 0.000001)

# ================================================================================== #

# Minimum function to determine when  the allele even get to the northern most populations
min_func <- function(x) {
if (all(is.na(x))) (NA) else min(x, na.rm = TRUE)
}

# Make scientific function more presentable
fancy_scientific <- function(l) {
  l <- as.character(l)
  l <- gsub("e", " * 0^", l)
  parse(text = l)
}

# ================================================================================== #

# Load data and merge files for each m

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m_SSP119/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m_SSP119/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m_SSP119 <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m_SSP119/phclineAFs_future_freq.')) %>% 
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
future_sim_data_vary_m_SSP119$m[which(future_sim_data_vary_m_SSP119$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m_SSP119$m[which(future_sim_data_vary_m_SSP119$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m_SSP119$m <- as.numeric(future_sim_data_vary_m_SSP119$m)
# Make column with scientific notation
future_sim_data_vary_m_SSP119 <- future_sim_data_vary_m_SSP119 %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m_SSP119, file = "data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP119.RData")
#load("data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP119.RData")


# Join
future_sim_data_vary_m_SSP119 %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m_SSP119$Site <- factor(future_sim_data_vary_m_SSP119$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Make m a factor
future_sim_data_vary_m_SSP119$m <- factor(future_sim_data_vary_m_SSP119$m, ms)

# Subset
future_sim_data_vary_m_SSP119_sub5 <- future_sim_data_vary_m_SSP119 %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP119_sub_5pop.pdf", width = 5.5, height = 14)
ggplot(future_sim_data_vary_m_SSP119_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()


# Graph AF vs time - per sim
future_sim_data_vary_m_SSP119_sub5_msub <- future_sim_data_vary_m_SSP119_sub5[which(future_sim_data_vary_m_SSP119_sub5$m == 0.01 | future_sim_data_vary_m_SSP119_sub5$m == 0.0001  | future_sim_data_vary_m_SSP119_sub5$m == 0.000001),]
# Average across the iterations
futureSSP119_avg_sub5_msub <- future_sim_data_vary_m_SSP119_sub5_msub %>% group_by(m, m_scientific, year, Site) %>% 
        reframe(AF_fut_avg = mean(AF_fut), AF_fut_sd = sd(AF_fut), n = n())
# Join
futureSSP119_avg_sub5_msub <- left_join(futureSSP119_avg_sub5_msub, ecovars, by = "Site")

# Graph subset
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP119_sub_5pop_msub.pdf", width = 5, height = 7)
ggplot() + 
        geom_line(data = future_sim_data_vary_m_SSP119_sub5_msub, aes(x = year+2000, y = AF_fut, color = Site, group=interaction(m, repId, Site)), 
                linewidth = 0.5, alpha = 0.01) + 
        geom_line(data = futureSSP119_avg_sub5_msub, aes(x = year+2000, y = AF_fut_avg, group = Site, color = Site), linewidth = 2) +
        facet_wrap(~m_scientific, ncol = 1, labeller = labeller(m_scientific = function(x) paste0("m = ", x))) + 
        scale_x_continuous(breaks = c(2020, 2060, 2100)) +
        scale_y_continuous(breaks = c(0, 0.5, 1)) +
        scale_color_manual(values = mycolors_sub5) + 
        labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + 
        theme(legend.position = "none",
            strip.text = element_text(size = 14, colour = "black", margin = margin(t = 3, b = 3, l = 0, r = 0)),
            strip.background = element_rect(fill = "white", colour = "black"), 
            plot.margin = margin(t = 40, r = 30,  b = 5.5, l = 5.5))
dev.off()

# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north_SSP119 <-
    foreach(mig=levels(future_sim_data_vary_m_SSP119$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m_SSP119 %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000,
        FC_common = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0.05)])+2000,
        avg_AF_2100 = mean(tmp2$AF_fut[which(tmp2$Site == "FC" & tmp2$year == 100)]))
    }
}
# Allele makes it to FC at high migration, but AF is really low

# Make numeric
future_sim_north_SSP119$m <- as.numeric(future_sim_north_SSP119$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_SSP119_summary <- future_sim_north_SSP119 %>%
  group_by(m) %>%
  summarise(arrived_count = sum(!is.na(FC)), arrived_perc = arrived_count/length(FC)*100, 
            common_count = sum(!is.na(FC_common)), common_perc = common_count/length(FC_common)*100)

# Make column with scientific notation
future_sim_north_SSP119_summary %<>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP119.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP119_summary, aes(y = na_perc, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations\n Allele Failed\n to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph - For each migration, how many sims failed to have allele become common at FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP119_common.pdf", width = 5.2714, height = 3.1)
ggplot(future_sim_north_SSP119_summary, aes(y = na_perc_common, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Sims\n MAF > 0.05") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none", plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
dev.off()

# Melt
future_sim_north_SSP119_summary_melt <- future_sim_north_SSP119_summary %>% 
    melt(id.vars = c("m", "m_scientific"), measure.vars = c("arrived_perc", "common_perc"), variable.name = "type", value.name = "perc")

# Graph - For each migration, how many sims failed to have allele get to FC or become common?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP119_arrived_common.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP119_summary_melt, aes(y = perc, x = m_scientific, shape = type)) + 
    geom_point(position = position_dodge(width = 0.5), size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m_SSP126/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m_SSP126/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m_SSP126 <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m_SSP126/phclineAFs_future_freq.')) %>% 
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
future_sim_data_vary_m_SSP126$m[which(future_sim_data_vary_m_SSP126$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m_SSP126$m[which(future_sim_data_vary_m_SSP126$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m_SSP126$m <- as.numeric(future_sim_data_vary_m_SSP126$m)
# Make column with scientific notation
future_sim_data_vary_m_SSP126 <- future_sim_data_vary_m_SSP126 %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m_SSP126, file = "data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP126.RData")
#load("data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP126.RData")


# Join
future_sim_data_vary_m_SSP126 %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m_SSP126$Site <- factor(future_sim_data_vary_m_SSP126$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Make m a factor
future_sim_data_vary_m_SSP126$m <- factor(future_sim_data_vary_m_SSP126$m, ms)

# Subset
future_sim_data_vary_m_SSP126_sub5 <- future_sim_data_vary_m_SSP126 %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP126_sub_5pop.pdf", width = 5.5, height = 14)
ggplot(future_sim_data_vary_m_SSP126_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) +  scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()


# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north_SSP126 <-
    foreach(mig=levels(future_sim_data_vary_m_SSP126$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m_SSP126 %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000,
        FC_common = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0.05)])+2000,
        avg_AF_2100 = mean(tmp2$AF_fut[which(tmp2$Site == "FC" & tmp2$year == 100)]))
    }
}

# Make numeric
future_sim_north_SSP126$m <- as.numeric(future_sim_north_SSP126$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_SSP126_summary <- future_sim_north_SSP126 %>%
  group_by(m) %>%
  summarise(na_count = sum(is.na(FC)), na_perc = na_count/length(FC)*100, 
            na_count_common = sum(is.na(FC_common)), na_perc_common = na_count_common/length(FC_common)*100)

# Make column with scientific notation
future_sim_north_SSP126_summary %<>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))


# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP126.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP126_summary, aes(y = na_perc, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations\n Allele Failed\n to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph - For each migration, how many sims failed to have allele become common at FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP126_common.pdf", width = 5.2714, height = 3.1)
ggplot(future_sim_north_SSP126_summary, aes(y = na_perc_common, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Sims\n MAF > 0.05") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none", plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
dev.off()

# ================================================================================== #

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m_SSP245/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m_SSP245/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m_SSP245 <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m_SSP245/phclineAFs_future_freq.')) %>% 
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
future_sim_data_vary_m_SSP245$m[which(future_sim_data_vary_m_SSP245$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m_SSP245$m[which(future_sim_data_vary_m_SSP245$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m_SSP245$m <- as.numeric(future_sim_data_vary_m_SSP245$m)
# Make column with scientific notation
future_sim_data_vary_m_SSP245 <- future_sim_data_vary_m_SSP245 %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m_SSP245, file = "data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP245.RData")
#load("data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP245.RData")


# Join
future_sim_data_vary_m_SSP245 %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m_SSP245$Site <- factor(future_sim_data_vary_m_SSP245$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Make m a factor
future_sim_data_vary_m_SSP245$m <- factor(future_sim_data_vary_m_SSP245$m, ms)

# Subset
future_sim_data_vary_m_SSP245_sub5 <- future_sim_data_vary_m_SSP245 %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP245_sub_5pop.pdf", width = 5.5, height = 14)
ggplot(future_sim_data_vary_m_SSP245_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()


# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north_SSP245 <-
    foreach(mig=levels(future_sim_data_vary_m_SSP245$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m_SSP245 %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000,
        FC_common = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0.05)])+2000,
        avg_AF_2100 = mean(tmp2$AF_fut[which(tmp2$Site == "FC" & tmp2$year == 100)]))
    }
}

# Make numeric
future_sim_north_SSP245$m <- as.numeric(future_sim_north_SSP245$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_SSP245_summary <- future_sim_north_SSP245 %>%
  group_by(m) %>%
  summarise(na_count = sum(is.na(FC)), na_perc = na_count/length(FC)*100, 
            na_count_common = sum(is.na(FC_common)), na_perc_common = na_count_common/length(FC_common)*100)

# Make column with scientific notation
future_sim_north_SSP245_summary %<>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))


# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP245.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP245_summary, aes(y = na_perc, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations\n Allele Failed\n to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph - For each migration, how many sims failed to have allele become common at FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP245_common.pdf", width = 5.2714, height = 3.1)
ggplot(future_sim_north_SSP245_summary, aes(y = na_perc_common, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Sims\n MAF > 0.05") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none", plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
dev.off()

# ================================================================================== #

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m_SSP370/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m_SSP370/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m_SSP370 <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m_SSP370/phclineAFs_future_freq.')) %>% 
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
future_sim_data_vary_m_SSP370$m[which(future_sim_data_vary_m_SSP370$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m_SSP370$m[which(future_sim_data_vary_m_SSP370$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m_SSP370$m <- as.numeric(future_sim_data_vary_m_SSP370$m)
# Make column with scientific notation
future_sim_data_vary_m_SSP370 <- future_sim_data_vary_m_SSP370 %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m_SSP370, file = "data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP370.RData")
#load("data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP370.RData")


# Join
future_sim_data_vary_m_SSP370 %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m_SSP370$Site <- factor(future_sim_data_vary_m_SSP370$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Make m a factor
future_sim_data_vary_m_SSP370$m <- factor(future_sim_data_vary_m_SSP370$m, ms)

# Subset
future_sim_data_vary_m_SSP370_sub5 <- future_sim_data_vary_m_SSP370 %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP370_sub_5pop.pdf", width = 5.5, height = 14)
ggplot(future_sim_data_vary_m_SSP370_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north_SSP370 <-
    foreach(mig=levels(future_sim_data_vary_m_SSP370$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m_SSP370 %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000,
        FC_common = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0.05)])+2000,
        avg_AF_2100 = mean(tmp2$AF_fut[which(tmp2$Site == "FC" & tmp2$year == 100)]))
    }
}

# Make numeric
future_sim_north_SSP370$m <- as.numeric(future_sim_north_SSP370$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_SSP370_summary <- future_sim_north_SSP370 %>%
  group_by(m) %>%
  summarise(na_count = sum(is.na(FC)), na_perc = na_count/length(FC)*100, 
            na_count_common = sum(is.na(FC_common)), na_perc_common = na_count_common/length(FC_common)*100)

# Make column with scientific notation
future_sim_north_SSP370_summary %<>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))


# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP370.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP370_summary, aes(y = na_perc, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations\n Allele Failed\n to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph - For each migration, how many sims failed to have allele become common at FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP370_common.pdf", width = 5.2714, height = 3.1)
ggplot(future_sim_north_SSP370_summary, aes(y = na_perc_common, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Sims\n MAF > 0.05") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none", plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
dev.off()

# ================================================================================== #

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_future/ph_vary_m_SSP585/', pattern = "phclineAFs_future_freq.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_future/ph_vary_m_SSP585/', x))))

# Loop through each N and extract files and merge
future_sim_data_vary_m_SSP585 <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  

      # State which file loading
      message(i)

      # Load file and add file name identifier
      tmp <- read.table(i) %>%
            mutate(file_name = i) %>%
            mutate(file_name = str_remove(file_name, pattern = 'data/processed/SLiM/ph_future/ph_vary_m_SSP585/phclineAFs_future_freq.')) %>% 
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
future_sim_data_vary_m_SSP585$m[which(future_sim_data_vary_m_SSP585$m == '1.0e-05')] <- "0.00001"
future_sim_data_vary_m_SSP585$m[which(future_sim_data_vary_m_SSP585$m == '1.0e-06')] <- "0.000001"
# Make column numeric
future_sim_data_vary_m_SSP585$m <- as.numeric(future_sim_data_vary_m_SSP585$m)
# Make column with scientific notation
future_sim_data_vary_m_SSP585 <- future_sim_data_vary_m_SSP585 %>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))

# Save output
save(future_sim_data_vary_m_SSP585, file = "data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP585.RData")
#load("data/processed/SLiM/ph_future/sim_future_data_vary_m_SSP585.RData")

# Join
future_sim_data_vary_m_SSP585 %<>% left_join(ecovars, by = "sim_eq")
# Make Site factor
future_sim_data_vary_m_SSP585$Site <- factor(future_sim_data_vary_m_SSP585$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Make m a factor
future_sim_data_vary_m_SSP585$m <- factor(future_sim_data_vary_m_SSP585$m, ms)

# Subset
future_sim_data_vary_m_SSP585_sub5 <- future_sim_data_vary_m_SSP585 %>% filter(Site == "FC" | Site == "SH" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph AF vs time - per sim
pdf("output/figures/SLiM/ph_future/AF_time_vary_m_SSP585_sub_5pop.pdf", width = 5.5, height = 14)
ggplot(future_sim_data_vary_m_SSP585_sub5, aes(x = (year+2000), y = AF_fut, color = Site, group=interaction(repId, Site))) + 
    geom_line(linewidth = 1, alpha = 0.2) + scale_y_continuous(breaks = c(0, 0.5, 1)) +
    facet_wrap(~m, ncol = 1) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Years", y = "Allele Frequency") + theme_linedraw(base_size = 30) + theme(strip.background = element_rect(fill = "white",colour = NA), strip.text = element_text(face="bold", color = "black"))+ theme(legend.position = "none")
dev.off()

# Loop through the iterations and calc when allele gets to northernmost pop
future_sim_north_SSP585 <-
    foreach(mig=levels(future_sim_data_vary_m_SSP585$m), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter for m
        tmp <- future_sim_data_vary_m_SSP585 %>% filter(m == mig)
        
        # Loop through each rep
        foreach(i=unique(tmp$repId), .combine="rbind", .errorhandling = "remove")%do%{  
        # Filter per rep
        tmp2 <- tmp %>% filter(repId == i)

        # Make table that calculates when northernmost pop first get the allele
        o =
        data.frame(
        repId = i,
        m = mig,
        FC = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0)])+2000,
        FC_common = min_func(tmp2$year[which(tmp2$Site == "FC" & tmp2$AF_fut > 0.05)])+2000,
        avg_AF_2100 = mean(tmp2$AF_fut[which(tmp2$Site == "FC" & tmp2$year == 100)]))
    }
}

# Make numeric
future_sim_north_SSP585$m <- as.numeric(future_sim_north_SSP585$m)

# How many sims didn't get the allele by 2100 in the northernmost pop?
future_sim_north_SSP585_summary <- future_sim_north_SSP585 %>%
  group_by(m) %>%
  summarise(na_count = sum(is.na(FC)), na_perc = na_count/length(FC)*100, 
            na_count_common = sum(is.na(FC_common)), na_perc_common = na_count_common/length(FC_common)*100)

# Make column with scientific notation
future_sim_north_SSP585_summary %<>% mutate(m_scientific = format(m, scientific = TRUE, digits = 1))


# Graph - For each migration, how many sims failed to have allele get to FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP585.pdf", width = 6, height = 3)
ggplot(future_sim_north_SSP585_summary, aes(y = na_perc, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Simulations\n Allele Failed\n to Arrive") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none")
dev.off()

# Graph - For each migration, how many sims failed to have allele become common at FC?
pdf("output/figures/SLiM/ph_future/AF_FC_pops_vary_m_SSP585_common.pdf", width = 5.2714, height = 3.1)
ggplot(future_sim_north_SSP585_summary, aes(y = na_perc_common, x = m_scientific)) + 
    geom_point(size = 5, color = mycolors_sub5[1]) +
    #scale_fill_manual(values =  mycolors_sub5[1]) + scale_color_manual(values =  mycolors_sub5[1]) + 
    scale_x_discrete(labels = fancy_scientific) + 
    scale_y_continuous(limits=c(-5, 105), breaks = c(0, 50, 100))+
    labs(x = "Migration Rate", y = "% Sims\n MAF > 0.05") + theme_linedraw(base_size = 30) + 
    guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
    theme(legend.position = "none", plot.margin = margin(t = 30, r = 5, b = 5, l = 5))
dev.off()
