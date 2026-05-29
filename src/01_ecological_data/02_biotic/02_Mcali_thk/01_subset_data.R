# Subset M. californianus thickness dataset

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'Hmisc', 'lmerTest'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(Hmisc)
library(lmerTest)

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/GEA/enviro_data/Mcali_thk")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_dir_fig <- paste("output/figures/enviro_data/Mcali_thk")
if (!dir.exists(out_dir_fig)) {dir.create(out_dir_fig)}

# ================================================================================== #

# Color palettes for graphing
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# ================================================================================== #

# Load data

# Read in mussel shell data
Mcali <- read.csv("data/raw/Mcali_thk/Mcali_shell_thk.csv", header=T)

# Read metadata
meta <- read.csv("guide_files/Populations_metadata.csv", header=T)
# Rename Site as Site.Code
meta <- rename(meta, Site.Code=Site)

# ================================================================================== #

# Create bin column
Mcali <- Mcali %>% mutate(
        value_bin = case_when(
          Length.L..mm. < 20 ~ "10",
          Length.L..mm. >= 20 & Length.L..mm. < 30 ~ "20",
          Length.L..mm. >= 30 & Length.L..mm. < 40 ~ "30",
          Length.L..mm. >= 40 & Length.L..mm. < 50 ~ "40",
          Length.L..mm. >= 50 & Length.L..mm. < 60 ~ "50",
          Length.L..mm. >= 60 & Length.L..mm. < 70 ~ "60",
          Length.L..mm. >= 70 & Length.L..mm. < 80 ~ "70",
          Length.L..mm. >= 80 & Length.L..mm. < 90 ~ "80",
          Length.L..mm. >= 90 & Length.L..mm. < 100 ~ "90",
          Length.L..mm. >= 100 & Length.L..mm. < 110 ~ "100",
          Length.L..mm. >= 110 & Length.L..mm. < 120 ~ "110",
          Length.L..mm. >= 120 & Length.L..mm. < 130 ~ "120",
          Length.L..mm. >= 130 & Length.L..mm. < 140 ~ "130",
          Length.L..mm. >= 140 & Length.L..mm. < 150 ~ "140",
          Length.L..mm. >= 150 ~ "150",
          TRUE ~ "Other" # Catch-all for values not fitting previous conditions
        )
      )

# ================================================================================== #

# Summarize for total number of mussel collected
Mcali.sum.collected <- Mcali %>% group_by(Site.Code, value_bin) %>% summarise(count=n())

# Write table
write.csv(Mcali.sum.collected, "data/raw/Mcali_thk/Mcali.sum.collected.csv", row.names=FALSE)

# ================================================================================== #

# Filter for only mussels that have done ImageJ on
Mcali_data <- Mcali %>% filter(!is.na(Segment.Area))

# Summarize for total number of mussel collected
Mcali.data.sum <- Mcali_data %>% group_by(Site.Code, value_bin) %>% summarise(count=n())
# Write table
write.csv(Mcali.data.sum, "data/raw/Mcali_thk/Mcali.data.sum.csv", row.names=FALSE)

# ================================================================================== #

# Set seed
set.seed(123)
#set.seed(179)

# Subsample dataset so even number of mussels of each size for each pop
Mcali_data_40 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "40") %>% slice_sample(n=3)
Mcali_data_50 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "50") %>% slice_sample(n=4)
Mcali_data_60 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "60") %>% slice_sample(n=5)
Mcali_data_70 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "70") %>% slice_sample(n=6)
Mcali_data_80 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "80") %>% slice_sample(n=5)
Mcali_data_90 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "90") %>% slice_sample(n=4)
Mcali_data_100 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "100") %>% slice_sample(n=3)

Mcali_data_sub <- rbind(Mcali_data_40, Mcali_data_50, Mcali_data_60, Mcali_data_70, Mcali_data_80, Mcali_data_90, Mcali_data_100)

# Add one KH mussel that is 80.8mm becuase 70-80mm is one short
Mcali_data_sub <- rbind(Mcali_data_sub, Mcali[1097,])

# ================================================================================== #

# Save dataset
write.csv(Mcali_data_sub, "data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_subset.csv", row.names=F)
# Read data
Mcali_data_sub <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_subset.csv")


# ================================================================================== #

# Format data

# Make sites a factors
Mcali_data_sub$Site.Code <- factor(Mcali_data_sub$Site.Code, 
    levels = rev(c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")))

# Make Integrated Thk numeric
Mcali_data_sub$Integrated.Thk <- as.numeric(Mcali_data_sub$Integrated.Thk)
# Make Max Thk numeric
Mcali_data_sub$Max.thk <- as.numeric(Mcali_data_sub$Max.thk)
# Make Min Thk numeric
Mcali_data_sub$Min.thk <- as.numeric(Mcali_data_sub$Min.thk)

# ================================================================================== #

# Summarize data
Mcali_sub_sum <- Mcali_data_sub %>%
    group_by(Site.Code) %>%
    summarise(num_shells = n(),
    mean_STI = mean(STI), sd_STI = sd(STI), se_STI = sd_STI/sqrt(num_shells), 
    mean_integrated_thk = mean(Integrated.Thk), sd_integrated_thk = sd(Integrated.Thk), se_integrated_thk = sd_integrated_thk/sqrt(num_shells), 
    mean_max_thk = mean(Max.thk), sd_max_thk = sd(Max.thk), se_max_thk = sd_max_thk/sqrt(num_shells),
    mean_min_thk = mean(Min.thk), sd_min_thk = sd(Min.thk), se_min_thk = sd_min_thk/sqrt(num_shells)) %>% as.data.frame() 

# ================================================================================== #

# Graph STI

# Graph by population - mean and se
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_STI_mean_se.pdf", width = 10, height = 14)
ggplot(data = Mcali_sub_sum, aes(x=mean_STI, y=Site.Code, colour = Site.Code)) + 
geom_point(data = Mcali_sub_sum, size=4, shape = 21, color = mycolors, fill = mycolors) + 
geom_errorbar(data = Mcali_sub_sum, aes(xmin=mean_STI-se_STI, xmax=mean_STI+se_STI), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlab("M. californianus Shell Thickness Index") + ylab("") + 
theme_classic(base_size = 24)
dev.off()

# Graph by population - raw points and pointrange
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_STI_raw_pointrange.pdf", width = 10, height = 14)
ggplot(data = Mcali_data_sub, aes(x=STI, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub, aes(x=STI, y=Site.Code), colour="darkgrey", height = 0.1) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color=mycolors) + 
theme_classic(base_size = 24) + xlim(0,6.5) +
xlab("M. californianus Shell Thickness Index") + ylab("")
dev.off()

# ================================================================================== #


# Graph Integrated thk

# Graph by population - mean and se
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_integrated_thk_mean_se.pdf", width = 10, height = 14)
ggplot(data = Mcali_sub_sum, aes(x=mean_integrated_thk, y=Site.Code, colour = Site.Code)) + 
geom_point(data = Mcali_sub_sum, size=4, shape = 21, color = mycolors, fill = mycolors) + 
geom_errorbar(data = Mcali_sub_sum, aes(xmin=mean_integrated_thk-se_integrated_thk, xmax=mean_integrated_thk+se_integrated_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlab("M. californianus Integrated Thickness (mm)") + ylab("") + 
theme_classic(base_size = 24)
dev.off()

# Graph by population - raw points and pointrange
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_integrated_thk_raw_pointrange.pdf", width = 10, height = 14)
ggplot(data = Mcali_data_sub, aes(x=Integrated.Thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub, aes(x=Integrated.Thk, y=Site.Code), colour="darkgrey", height = 0.1) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color=mycolors) + 
theme_classic(base_size = 24) + xlim(0, 5) +
xlab("M. californianus Integrated Thickness (mm)") + ylab("")
dev.off()

# ================================================================================== #

# Graph Max thk

# Graph by population - mean and se
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_max_thk_mean_se.pdf", width = 10, height = 14)
ggplot(data = Mcali_sub_sum, aes(x=mean_max_thk, y=Site.Code, colour = Site.Code)) + 
geom_point(data = Mcali_sub_sum, size=4, shape = 21, color = mycolors, fill = mycolors) + 
geom_errorbar(data = Mcali_sub_sum, aes(xmin=mean_max_thk-se_max_thk, xmax=mean_max_thk+se_max_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlab("M. californianus Max Thickness (mm)") + ylab("") + 
theme_classic(base_size = 24)
dev.off()

# Graph by population - raw points and pointrange
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_max_thk_raw_pointrange.pdf", width = 10, height = 14)
ggplot(data = Mcali_data_sub, aes(x=Max.thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub, aes(x=Max.thk, y=Site.Code), colour="darkgrey", height = 0.1) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color=mycolors) + 
theme_classic(base_size = 24) + xlim(0,7) +
xlab("M. californianus Max Thickness (mm)") + ylab("")
dev.off()

# ================================================================================== #

# Graph Min thk

# Graph by population - mean and se
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_min_thk_mean_se.pdf", width = 10, height = 14)
ggplot(data = Mcali_sub_sum, aes(x=mean_min_thk, y=Site.Code, colour = Site.Code)) + 
geom_point(data = Mcali_sub_sum, size=4, shape = 21, color = mycolors, fill = mycolors) + 
geom_errorbar(data = Mcali_sub_sum, aes(xmin=mean_min_thk-se_min_thk, xmax=mean_min_thk+se_min_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlab("M. californianus Min Thickness (mm)") + ylab("") + 
theme_classic(base_size = 24)
dev.off()

# Graph by population - raw points and pointrange
pdf("output/figures/GEA/enviro/Mcali_thk/Mcali_min_thk_raw_pointrange.pdf", width = 10, height = 14)
ggplot(data = Mcali_data_sub, aes(x=Min.thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub, aes(x=Min.thk, y=Site.Code), colour="darkgrey", height = 0.1) +
stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color=mycolors) + 
theme_classic(base_size = 24) + xlim(0,3) +
xlab("M. californianus Min Thickness (mm)") + ylab("")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Format data for GLMs
Mcalifornianus_data <- left_join(meta, Mcali_sub_sum, by="Site.Code")

# For glms only use means of thickness metrics
Mcalifornianus_data_clean <- Mcalifornianus_data[, c(1:5,6,9,12,15)]

# ================================================================================== #

# Save dataset
write.csv(Mcalifornianus_data_clean, "data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean.csv", row.names=F)

# ================================================================================== #
# ================================================================================== #

# Graphing as map projections

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Order Metadata
Mcali_data_sub$Site.Code <- factor(Mcali_data_sub$Site.Code, 
levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD","KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Summarize data
Mcali_data_sub_sum <- Mcali_data_sub %>% group_by(Site.Code) %>% summarise(mean_STI = mean(STI), mean_integrated_thk = mean(Integrated.Thk))
# Bind
Mcali_data_sub_sum_meta <- left_join(Mcali_data_sub_sum, meta, by="Site.Code")

# Cut ARA since outlier for collection
Mcali_data_sub_sum.18 <- Mcali_data_sub_sum_meta[-which(Mcali_data_sub_sum_meta$Site.Code == "ARA"),]

# Graph STI 2024 data
pdf("output/figures/enviro/Mcali_thk/Mcali_STI_2024_18sites_altcolors.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_STI), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  #scale_fill_viridis(option="viridis", direction = -1) +
  scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name="STI") +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 24) + #ggtitle("STI Shell Thickness Projections") +
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = c(0.98, 0.52))
dev.off()

# Graph Mean Integrated Thk 2024 data
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_altcolors.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name="Integrated Thickness") +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 24) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 16), legend.position = c(0.98, 0.52))
dev.off()
# Graph Mean Integrated Thk 2024 data (alt)
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_altsize.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name=NULL) +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw(base_size = 26) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.80, 0.53), legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
dev.off()
# Graph Mean Integrated Thk 2024 data (alt colors)
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_altsizecol.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlGnBu"), name=NULL) +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw(base_size = 26) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.80, 0.53), legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
dev.off()
# Graph Mean Integrated Thk 2024 data (biger)
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_bigger.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlGnBu"), name="Thickness") +
             coord_fixed(1.3) +
  xlim(c(-125, -112)) +
  #ylim(c(32, 46)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw(base_size = 26) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.80, 0.47), legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
dev.off()