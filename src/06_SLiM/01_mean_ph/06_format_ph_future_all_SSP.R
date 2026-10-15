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
mycolors_sub <- hcl.colors(n = 7, palette = "SunsetDark")[c(1,3,4,6)]

# ================================================================================== #

# Generate output directories

# Bio-Oracle directory
out_dir_future <- paste("data/raw/Bio-oracle/future")
if (!dir.exists(out_dir_future)) {dir.create(out_dir_future)}

# Data SliM directory
out_data_dir <- paste("data/processed/SLiM/ph_future")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# Load data
bio_oracle_ssp585 <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp585.csv")
write.csv(bio_oracle_ssp585, "data/processed/enviro_data/Bio-oracle/bio_oracle_ssp585.csv", row.names=FALSE)

bio_oracle_ssp585_sites <- read.csv("data/processed/enviro_data/Bio-oracle/bio_oracle_ssp585_sites.csv")
write.csv(bio_oracle_ssp585_sites, "data/processed/enviro_data/Bio-oracle/bio_oracle_ssp585_sites.csv", row.names=FALSE)

bio_oracle_ssp585_sites

# Extract just mean pH


# Add decade term
ph_ssp585_sites %<>% group_by(location) %>% mutate(decade=seq(5,75, by=10))
 
# Calculate slope and intercept for each location
ph_future_lm <- foreach(i=unique(ph_ssp585_sites$location), .combine="rbind", .errorhandling = "remove")%do%{  
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
ph_future_lm$Site <- factor(ph_future_lm$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Order df by site 
ph_future_lm <- ph_future_lm[order(ph_future_lm$Site),]

# Save slope and intercept for each pop
write.csv(ph_future_lm, "data/processed/SLiM/ph_future/ph_future_lm.csv", row.names=F)

# ================================================================================== #

# Extract last decade 
ph_ssp585_sites_2090_tmp <- ph_ssp585_sites[which(ph_ssp585_sites$decade == 75),]
# Reverse order
ph_ssp585_sites_2090 <- ph_ssp585_sites_2090_tmp[nrow(ph_ssp585_sites_2090_tmp):1,]

write.csv(ph_ssp585_sites_2090, "data/processed/SLiM/ph_future/ph_future_2090_sites.csv", row.names=F)

# ================================================================================== #

# For visualization purposes make x-axis go from 2020-2100

# Calc pred pH adjusted for 2020-2100
pred <- ph_future_lm %>%rowwise() %>%
do(data.frame(Site = .$Site, year = 2020:2100, ph = .$intercept + .$slope * (2020:2100 - 2020)))

# Threshold from simulations
thresh = 7.991

# Graph
pdf("output/figures/SLiM/ph_future/pH_future_Bio_oracle.pdf", width = 10, height = 8)
ggplot(ph_ssp585_sites, aes(x = decade+2020, y = ph_mean, color = location)) + geom_point(size = 5, shape = 16, alpha = 0.9) +
    geom_line(data = pred, aes(x = year, y = ph, color = Site), linewidth = 1)+
    #geom_abline(data = ph_future_lm, aes(slope = slope, intercept = intercept, color = Site)) + 
    geom_hline(yintercept=thresh, col = "black", linetype = "dashed") +
    ylim(7.6, 8.078) +
    scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + 
    labs(x = "Year", y = "Mean pH") + theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()
