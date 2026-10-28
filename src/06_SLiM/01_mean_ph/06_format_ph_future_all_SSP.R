# Format pH future Bio-oracle data

# Clear memory
rm(list=ls())

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'devtools', 'magrittr, 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(tidyverse)
library(devtools) 
#devtools::install_github("bio-oracle/biooracler")
library(biooracler)
library(magrittr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub5 <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6,7)]

# ================================================================================== #

# Generate output directories

# Bio-Oracle directory
out_dir_future <- paste("data/raw/Bio-oracle/future")
if (!dir.exists(out_dir_future)) {dir.create(out_dir_future)}

# Data SliM directory
out_data_dir <- paste("data/processed/SLiM/ph_future")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# Threshold from simulations
thresh = 7.991

# ================================================================================== #

# Load data
bio_oracle_ssp585_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp585_sites.csv")

# Extract just mean pH
ph_ssp585_sites <- bio_oracle_ssp585_sites[,c(1,2,3,11,13)]

# Add decade term
ph_ssp585_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_ssp585_lm <- foreach(i=unique(ph_ssp585_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp585_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp585_sites$location <- factor(ph_ssp585_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_ssp585_lm$Site <- factor(ph_future_ssp585_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_ssp585_lm <- ph_future_ssp585_lm[order(ph_future_ssp585_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_ssp585_lm, "data/processed/SLiM/ph_future/ph_future_ssp585_lm.csv", row.names=F)


# For visualization purposes make x-axis go from 2020-2100
# Calc pred pH adjusted for 2020-2100
pred_ssp585 <- ph_future_ssp585_lm %>%rowwise() %>% do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp585.pdf", width = 5.5, height = 5.5)
ggplot(ph_ssp585_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp585, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph just 5 simulated sites
ph_ssp585_sites5 <- ph_ssp585_sites %>% filter(location == "FC" | location == "ARA" | location == "CBL" | location == "BMR" | location == "STR")
pred_ssp585sites5 <- pred_ssp585 %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp585_sub_5pop.pdf", width = 5.98, height = 3)
ggplot(ph_ssp585_sites5, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp585sites5, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    scale_x_continuous(breaks = c(2020, 2060, 2100)) +
    scale_y_continuous(limits = c(7.6, 8.078), breaks = c(7.6, 7.8, 8.0)) +
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    scale_fill_manual(values = mycolors_sub5) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none") +
    theme(plot.margin = margin(t = 5.5, r = 30,  b = 5.5, l = 5.5))
dev.off()

# ================================================================================== #

# Load data
bio_oracle_ssp370_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp370_sites.csv")

# Extract just mean pH
ph_ssp370_sites <- bio_oracle_ssp370_sites[,c(1,2,3,11,13)]

# Add decade term
ph_ssp370_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_ssp370_lm <- foreach(i=unique(ph_ssp370_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp370_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp370_sites$location <- factor(ph_ssp370_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_ssp370_lm$Site <- factor(ph_future_ssp370_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_ssp370_lm <- ph_future_ssp370_lm[order(ph_future_ssp370_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_ssp370_lm, "data/processed/SLiM/ph_future/ph_future_ssp370_lm.csv", row.names=F)


# For visualization purposes make x-axis go from 2020-2100
# Calc pred pH adjusted for 2020-2100
pred_ssp370 <- ph_future_ssp370_lm %>%rowwise() %>% do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp370.pdf", width = 5.5, height = 5.5)
ggplot(ph_ssp370_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp370, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()


# ================================================================================== #

# Load data
bio_oracle_ssp245_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp245_sites.csv")

# Extract just mean pH
ph_ssp245_sites <- bio_oracle_ssp245_sites[,c(1,2,3,11,13)]

# Add decade term
ph_ssp245_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_ssp245_lm <- foreach(i=unique(ph_ssp245_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp245_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp245_sites$location <- factor(ph_ssp245_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_ssp245_lm$Site <- factor(ph_future_ssp245_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_ssp245_lm <- ph_future_ssp245_lm[order(ph_future_ssp245_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_ssp245_lm, "data/processed/SLiM/ph_future/ph_future_ssp245_lm.csv", row.names=F)


# For visualization purposes make x-axis go from 2020-2100
# Calc pred pH adjusted for 2020-2100
pred_ssp245 <- ph_future_ssp245_lm %>%rowwise() %>% do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp245.pdf", width = 5.5, height = 5.5)
ggplot(ph_ssp245_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp245, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph just 5 simulated sites
ph_ssp245_sites5 <- ph_ssp245_sites %>% filter(location == "FC" | location == "ARA" | location == "CBL" | location == "BMR" | location == "STR")
pred_ssp245sites5 <- pred_ssp245 %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp245_sub_5pop.pdf", width = 5.98, height = 3)
ggplot(ph_ssp245_sites5, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp245sites5, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    scale_x_continuous(breaks = c(2020, 2060, 2100)) +
    scale_y_continuous(limits = c(7.6, 8.078), breaks = c(7.6, 7.8, 8.0)) +
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    scale_fill_manual(values = mycolors_sub5) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none") +
    theme(plot.margin = margin(t = 5.5, r = 30,  b = 5.5, l = 5.5))
dev.off()

# ================================================================================== #

# Load data
bio_oracle_ssp126_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp126_sites.csv")

# Extract just mean pH
ph_ssp126_sites <- bio_oracle_ssp126_sites[,c(1,2,3,11,13)]

# Add decade term
ph_ssp126_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_ssp126_lm <- foreach(i=unique(ph_ssp126_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp126_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp126_sites$location <- factor(ph_ssp126_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_ssp126_lm$Site <- factor(ph_future_ssp126_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_ssp126_lm <- ph_future_ssp126_lm[order(ph_future_ssp126_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_ssp126_lm, "data/processed/SLiM/ph_future/ph_future_ssp126_lm.csv", row.names=F)


# For visualization purposes make x-axis go from 2020-2100
# Calc pred pH adjusted for 2020-2100
pred_ssp126 <- ph_future_ssp126_lm %>%rowwise() %>% do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp126.pdf", width = 5.5, height = 5.5)
ggplot(ph_ssp126_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp126, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()


# ================================================================================== #

# Load data
bio_oracle_ssp119_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp119_sites.csv")

# Extract just mean pH
ph_ssp119_sites <- bio_oracle_ssp119_sites[,c(1,2,3,11,13)]

# Add decade term
ph_ssp119_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_ssp119_lm <- foreach(i=unique(ph_ssp119_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- ph_ssp119_sites %>% filter(location == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(ph_mean ~ decade, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2])
}
# Set site to factor
ph_ssp119_sites$location <- factor(ph_ssp119_sites$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
ph_future_ssp119_lm$Site <- factor(ph_future_ssp119_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_ssp119_lm <- ph_future_ssp119_lm[order(ph_future_ssp119_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_ssp119_lm, "data/processed/SLiM/ph_future/ph_future_ssp119_lm.csv", row.names=F)


# For visualization purposes make x-axis go from 2020-2100
# Calc pred pH adjusted for 2020-2100
pred_ssp119 <- ph_future_ssp119_lm %>%rowwise() %>% do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp119.pdf", width = 5.5, height = 5.5)
ggplot(ph_ssp119_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp119, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph just 5 simulated sites
ph_ssp119_sites5 <- ph_ssp119_sites %>% filter(location == "FC" | location == "ARA" | location == "CBL" | location == "BMR" | location == "STR")
pred_ssp119sites5 <- pred_ssp119 %>% filter(Site == "FC" | Site == "ARA" | Site == "CBL" | Site == "BMR" | Site == "STR")
# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle_ssp119_sub_5pop.pdf", width = 5.98, height = 3)
ggplot(ph_ssp119_sites5, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred_ssp119sites5, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    scale_x_continuous(breaks = c(2020, 2060, 2100)) +
    scale_y_continuous(limits = c(7.6, 8.078), breaks = c(7.6, 7.8, 8.0)) +
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    scale_fill_manual(values = mycolors_sub5) + scale_color_manual(values = mycolors_sub5) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none") +
    theme(plot.margin = margin(t = 5.5, r = 30,  b = 5.5, l = 5.5))
dev.off()
