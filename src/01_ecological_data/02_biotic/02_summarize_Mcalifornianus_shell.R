# Subset and summarize M. californianus thickness dataset

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set seed for reproducibility
set.seed(123)

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'Hmisc', 'lmerTest'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
viridiscolors <- viridis(n=19)
viridiscolors18 <- viridiscolors[-4]

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/GEA/enviro/Mcali_thk")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_dir_fig <- paste("output/figures/enviro/Mcali_thk")
if (!dir.exists(out_dir_fig)) {dir.create(out_dir_fig)}

# ================================================================================== #

# Load data

# Read in mussel shell data
Mcali <- read.csv("data/raw/Mcali_thk/Mcalifornianus_shell_traits.csv", header=T)

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

# Subsample dataset so even number of mussels of each size for each pop 
# Note: make sure set seed above for reproducibility
Mcali_data_40 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "40") %>% slice_sample(n=3)
Mcali_data_50 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "50") %>% slice_sample(n=4)
Mcali_data_60 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "60") %>% slice_sample(n=5)
Mcali_data_70 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "70") %>% slice_sample(n=6)
Mcali_data_80 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "80") %>% slice_sample(n=5)
Mcali_data_90 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "90") %>% slice_sample(n=4)
Mcali_data_100 <- Mcali_data %>% group_by(Site.Code) %>% filter(value_bin == "100") %>% slice_sample(n=3)

# Join
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

# Remove ARA since mussels collected at this site were in a more wave protected area
Mcali_data_sub_noARA <- Mcali_data_sub %>% filter(Site.Code != "ARA")
Mcali_sub_sum_noARA <- Mcali_sub_sum %>% filter(Site.Code != "ARA")

Mcali_data_sub_noARA_thk_sum <- Mcali_data_sub_noARA %>% 
  group_by(Site.Code) %>%
    summarise(mean_integrated_thk = mean(Integrated.Thk), max_integrated_thk = max(Integrated.Thk), min_integrated_thk = min(Integrated.Thk))

# Add column that highlights N, and S
Mcali_data_sub_noARA_thk_sum <- Mcali_data_sub_noARA_thk_sum %>% mutate(shape = case_when(Site.Code %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   Site.Code %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))

# Graph
pdf("output/figures/enviro/Mcali_thk/Mcali_thk_sites.pdf", width = 6, height = 7.18)
ggplot(Mcali_data_sub_noARA_thk_sum, aes(x = mean_integrated_thk, y = Site.Code, shape = shape)) +
  #geom_jitter(data = Mcali_data_sub_noARA, aes(x=Integrated.Thk, y=Site.Code), colour="darkgrey", height = 0.1, size = 2, alpha = 0.8) +
  geom_linerange(data = Mcali_data_sub_noARA_thk_sum, aes(xmin = min_integrated_thk, xmax = max_integrated_thk)) +
  scale_shape_manual(values = c(21, 23)) +
  geom_point(size = 6, color = "black", fill = rev(viridiscolors18)) +
  labs(x = "Thickness", y = "") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Graph raw points, means and sd

# Graph STI
pdf("output/figures/enviro/Mcali_thk/Mcali_STI_mean_sd.pdf", width = 7, height = 10)
ggplot(data = Mcali_sub_sum_noARA, aes(x=mean_STI, y=Site.Code,)) + 
geom_jitter(data = Mcali_data_sub_noARA, aes(x=STI, y=Site.Code), colour="darkgrey", height = 0.1, size = 3, alpha = 0.8) +
geom_point(data = Mcali_sub_sum_noARA, size=6) + 
geom_errorbar(data = Mcali_sub_sum_noARA, aes(xmin=mean_STI-sd_STI, xmax=mean_STI+sd_STI), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlim(0,7) +
xlab("STI") + ylab("") + 
theme_classic(base_size = 30)
dev.off()

# Graph Integrated thk
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thk_mean_sd.pdf", width = 7, height = 10)
ggplot(data = Mcali_sub_sum_noARA, aes(x=mean_integrated_thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub_noARA, aes(x=Integrated.Thk, y=Site.Code), colour="darkgrey", height = 0.1, size = 3, alpha = 0.8) +
geom_point(data = Mcali_sub_sum_noARA, size=6) + 
geom_errorbar(data = Mcali_sub_sum_noARA, aes(xmin=mean_integrated_thk-sd_integrated_thk, xmax=mean_integrated_thk+sd_integrated_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlim(0,7) +
xlab("Cross-sectional Thickness") + ylab("") + 
theme_classic(base_size = 30)
dev.off()

# Graph Max thk
pdf("output/figures/enviro/Mcali_thk/Mcali_max_thk_mean_sd.pdf", width = 7, height = 10)
ggplot(data = Mcali_sub_sum_noARA, aes(x=mean_max_thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub_noARA, aes(x=Max.thk, y=Site.Code), colour="darkgrey", height = 0.1, size = 3, alpha = 0.8) +
geom_point(data = Mcali_sub_sum_noARA, size=6) + 
geom_errorbar(data = Mcali_sub_sum_noARA, aes(xmin=mean_max_thk-sd_max_thk, xmax=mean_max_thk+sd_max_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlim(0,7) +
xlab("Max Thickness") + ylab("") + 
theme_classic(base_size = 30)
dev.off()

# Graph Min thk
pdf("output/figures/enviro/Mcali_thk/Mcali_min_thk_mean_sd.pdf", width = 7, height = 10)
ggplot(data = Mcali_sub_sum_noARA, aes(x=mean_min_thk, y=Site.Code)) + 
geom_jitter(data = Mcali_data_sub_noARA, aes(x=Min.thk, y=Site.Code), colour="darkgrey", height = 0.1, size = 3, alpha = 0.8) +
geom_point(data = Mcali_sub_sum_noARA, size=6) + 
geom_errorbar(data = Mcali_sub_sum_noARA, aes(xmin=mean_min_thk-sd_min_thk, xmax=mean_min_thk+sd_min_thk), width=.5) +
scale_fill_manual(values=mycolors, guide="none") + 
scale_colour_manual(values=mycolors, guide="none") +
xlim(0,7) +
xlab("Min Thickness") + ylab("") + 
theme_classic(base_size = 30)
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
# Join
Mcali_data_sub_sum_meta <- left_join(Mcali_data_sub_sum, meta, by="Site.Code")

# Cut ARA since outlier for collection
Mcali_data_sub_sum.18 <- Mcali_data_sub_sum_meta[-which(Mcali_data_sub_sum_meta$Site.Code == "ARA"),]


# Graph Mean Integrated Thk 2024 data
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

# Graph Mean Integrated Thk 2024 data (bigger)
pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_bigger.pdf", width = 9, height = 9)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 12) + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlGnBu"), name="Thickness") +
             coord_fixed(1.3) +
  #xlim(c(-125, -112.5)) +
   scale_x_continuous(limits = c(-125, -113), breaks = seq(-125, -113, by = 3)) +
  #ylim(c(32, 46)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 28) + #ggtitle("Integrated Shell Thickness Projections") +
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + 
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 22), legend.position = c(0.78, 0.515))
dev.off()

pdf("output/figures/enviro/Mcali_thk/Mcali_integrated_thick_18sites_2024_wider.pdf", width = 10, height = 10.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = Mcali_data_sub_sum.18, aes(x = Long, y = Lat, fill = mean_integrated_thk), shape = 21, size = 15) + 
  #scale_fill_gradientn(colours=rev(brewer.pal(9, "BrBG")), name="Thickness", breaks = c(1.6, 1.9, 2.2)) +
  scale_fill_viridis(option="mako", name="Thickness", breaks = c(1.6, 1.9, 2.2), direction = -1) +
             coord_fixed(1.3) +
   scale_x_continuous(limits = c(-125, -113), breaks = seq(-125, -113, by = 5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + #ggtitle("Integrated Shell Thickness Projections") +
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + 
  theme(legend.title = element_text(size = 30), legend.text = element_text(size = 28), legend.position = c(0.82, 0.44), legend.margin = margin(0, 0, 0, 0))
dev.off()