# Format Mcali data

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
library(magrittr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_sub <- mycolors[c(1,3,9,11,14)]

# ================================================================================== #

# Generate output directories

# Data directory
out_data_dir <- paste("data/processed/SLiM/Mcali_future")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/Mcali_future")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

# Shell thickness data from 2001-2019 (Longman & Sanford 2025)
Mcal_thk <- read.csv("data/raw/Mcali_thk/Mussel_Shell_Thickness_2001_2019.csv")

# Add new M cali data
Mcali_new <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_subset.csv")

# ================================================================================== #

# Format and clean data - 2001-2019

# Format
Mcal_thk$Site.Code <- factor(Mcal_thk$Site.Code, levels = c("FC", "SH", "ARA", "VD", "BMR", "SBR"))
Mcal_thk %<>% rename(year = Date.1, Site.full = Site, Site = Site.Code)

# Remove broken shells
Mcal_thk <- Mcal_thk[-which(is.na(Mcal_thk$Average.Thickness)),]
# Remove outlier
Mcal_thk <- Mcal_thk[-which(Mcal_thk$Time == "2008/09" & Mcal_thk$Length<35),]

# Calculate Log of thickness/shell length ratio
Mcal_thk$thk.length <- Mcal_thk$Average.Thickness/Mcal_thk$Length
Mcal_thk$log.thk.length <- log(Mcal_thk$thk.length)

# Drop ARA
Mcal_thk %<>% filter(Site != "ARA")

# Rescale x-axis so 2020 is 0
Mcal_thk %<>% mutate(year.scaled = year-2020)

# Calculate a thickness metric that will be comparable to current data which is adjusted based on size
Mcal_thk %<>% mutate(Thickness = thk.length*75)

# Subset to only a few column to join with previous data
Mcal_thk_sub <- Mcal_thk[,c("Site", "year", "year.scaled", "Thickness")]

# ================================================================================== #

# Format and clean data - 2001-2019

# Separate collection date
Mcali_new %<>% separate_wider_delim(cols = Date.Collected, delim = ".", names = c("month", "day", "year"))

# Format
Mcali_new %<>% rename(Site = Site.Code, Length = 'Length.L..mm.')
Mcali_new$year <- as.numeric(Mcali_new$year)

# Filter so only 5 focal sites
Mcali_new_filter <- Mcali_new %>% filter(Site %in% unique(Mcal_thk$Site))

# Calc avg of max and min
Mcali_new_filter %<>% mutate(Average.Thickness = (Max.thk+Min.thk)/2)

# Rescale x-axis so 2020 is 0
Mcali_new_filter %<>% mutate(year.scaled = year-2020)

# Calculate a thickness metric that will be comparable to current data which is adjusted based on size
Mcali_new_filter %<>% mutate(Thickness = Average.Thickness/Length*75)

# Subset to only a few column to join with previous data
Mcali_new_filter_sub <- Mcali_new_filter[,c("Site", "year", "year.scaled", "Thickness")]

# ================================================================================== #

# Join data
Mcali_all <- rbind(Mcal_thk_sub, Mcali_new_filter_sub)

# ================================================================================== #

# Calculate slope and intercept for each location
Mcali_lm <- foreach(i=unique(Mcali_all$Site), .combine="rbind", .errorhandling = "remove")%do%{  
    # Filter for specific site
    tmp <- Mcali_all %>% filter(Site == i)
    # Linear model - reformat time so centered on each decade
    mod <- lm(Thickness ~ year.scaled, tmp)
    # Format data
    data.frame(
    Site = i,
    intercept = mod$coefficients[1],
    slope = mod$coefficients[2], 
    p = summary(mod)$coefficients[2,4])
}

# Set site to factor
Mcali_lm$Site <- factor(Mcali_lm$Site, levels=c("FC", "SH", "VD", "BMR", "SBR"))
# Order df by site 
Mcali_lm <- Mcali_lm[order(Mcali_lm$Site),]

# ================================================================================== #


# Graph
pdf("output/figures/SLiM/Mcali_future/Mcali_lm.pdf", width = 6.5, height = 6)
ggplot(Mcali_all, aes(x = year.scaled, y = Thickness, fill = Site, color = Site)) + geom_jitter(size = 1, shape = 21, alpha = 0.9) +
    geom_abline(data = Mcali_lm, aes(slope = slope, intercept = intercept, color = Site), linewidth = 2) + 
    scale_fill_manual(values = mycolors_sub) + scale_color_manual(values = mycolors_sub) + 
    labs(x = "Years", y = "Shell Thk") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/SLiM/Mcali_future/Mcali_lm_2100.pdf", width = 6.5, height = 6)
ggplot(Mcali_all, aes(x = year.scaled, y = Thickness, fill = Site, color = Site)) + geom_jitter(size = 1, shape = 21, alpha = 0.9) +
    geom_abline(data = Mcali_lm, aes(slope = slope, intercept = intercept, color = Site), linewidth = 2) + 
    scale_fill_manual(values = mycolors_sub) + scale_color_manual(values = mycolors_sub) + 
    scale_x_continuous(limits = c(-21, 80), breaks = c(-20, 0, 20, 40, 60, 80)) + ylim(0, 3.5) +
    labs(x = "Years", y = "Shell Thk") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()





# For all sites calc avg thickness to make sure very similar to the integrate thk

afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T)


Mcali_new_test <-  Mcali_new %>% mutate(Average.Thickness = (Max.thk+Min.thk)/2)

# Summarize
Mcali_new_test_sum <- Mcali_new_test %>% group_by(Site) %>%
    summarise(mean_avg_thk = mean(Average.Thickness))

Mcali_compare <- left_join(afs.Mcali, Mcali_new_test_sum)
Mcali_compare <- arrange(Mcali_compare, desc(latitude))