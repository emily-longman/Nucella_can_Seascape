# ------------------------------------------------------------------------------
# Mussel collections estimate:  
# Determine the number of mussels and size dist needed to get a representative sample
# Emily K. Longman
# ------------------------------------------------------------------------------

# Set path as main Github repo
# Install and load package
#install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(maps) 
library(mapdata)
library(RColorBrewer)

# Load data --------------------------------------------------------------------

# Read in metadata 
metadata <- read.csv("guide_files/Populations_metadata.csv", header=T)
# Filter metadata for 2019 sites
metadata_sub <- metadata %>% filter(Site %in% c("FC", "SH", "ARA", "VD", "BMR", "SBR"))

# Load 2019 mussel data
mussels <- read.csv("data/raw/Mcali_thk/Mcal_shell_thickness_2019.csv")
mussels$Thk.length <- mussels$Ave.thk.1.3/ mussels$Shell.Length
str(mussels)

# Load 2023/2024 data
mussels_2024 <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_subset.csv", header=T)

# Look at all data -------------------------------------------------------------

# Graph all mussels
ggplot(mussels, aes(Shell.Length, Ave.thk.1.3)) + geom_point()

#Linear regression
mussel.mod <- lm(Ave.thk.1.3 ~ Shell.Length, mussels)
plot(mussel.mod)
summary(mussel.mod)

#calc thickness/ length ratio
mean(mussels$Thk.length) #0.02117015

#Just Bodega Mussels -----------------------------------------------------------
mussels.BMR <- mussels[which(mussels$Site.Code == "BH"),]

length(mussels.BMR) #101 mussels

#calc thickness/ length ratio
mean(mussels.BMR$Thk.length) #0.02124721

ggplot(mussels.BMR, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.BMR <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.BMR)
summary(mussel.mod.BMR)

#Sample just 50 BMR mussels across all sizes 100 times 
mussels.BMR.sample = NULL
for (i in 1:100) {
  s = sample(mussels.BMR$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.BMR.sample<- c(mussels.BMR.sample, m)
  }
hist(mussels.BMR.sample)
abline(v = 0.02124721, col= 2)

range(mussels.BMR.sample)


#Sample 50 BMR mussels that are less than 120mm 100 times 
mussels.BMR.less.120 <- mussels.BMR[which(mussels.BMR$Shell.Length < 120),]

mussels.BMR.less.120.sample = NULL
for (i in 1:100) {
  s = sample(mussels.BMR.less.120$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.BMR.less.120.sample<- c(mussels.BMR.less.120.sample, m)
  }
hist(mussels.BMR.less.120.sample)
abline(v = 0.02124721, col= 2)

range(mussels.BMR.less.120.sample)


# Just Strawberry Hill Mussels -------------------------------------------------
mussels.SH <- mussels[which(mussels$Site.Code == "SH"),]

dim(mussels.SH) #108 mussels
mean(mussels.SH$Ave.thk.1.3/ mussels.SH$Shell.Length) #0.02396629

ggplot(mussels.SH, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.SH <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.SH)
plot(mussel.mod.SH)
summary(mussel.mod.SH)

#Sample just 50 SH mussels across all sizes 100 times 
mussels.SH.sample = NULL
for (i in 1:100) {
  s = sample(mussels.SH$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.SH.sample<- c(mussels.SH.sample, m)
  }
hist(mussels.SH.sample)
abline(v = 0.02396629, col= 2)

range(mussels.SH.sample) #0.02275464 0.02591909

#Sample 50 SH mussels that are less than 120mm and greater than 40mm 100 times 
mussels.SH.less.120 <- mussels.SH[which(mussels.SH$Shell.Length < 120 & mussels.SH$Shell.Length > 40),]
mussels.SH.less.120.sample = NULL
for (i in 1:100) {
  s = sample(mussels.SH.less.120$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.SH.less.120.sample<- c(mussels.SH.less.120.sample, m)
  }
hist(mussels.SH.less.120.sample)
abline(v = 0.02396629, col= 2)

range(mussels.SH.less.120.sample)

# Just Fogarty Creek Mussels ---------------------------------------------------
mussels.FC <- mussels[which(mussels$Site.Code == "FC"),]

dim(mussels.FC) #108 mussels
mean(mussels.FC$Ave.thk.1.3/ mussels.FC$Shell.Length) #0.01994755

ggplot(mussels.FC, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.FC <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.FC)
plot(mussel.mod.FC)
summary(mussel.mod.FC)

#Sample just 50 SH mussels across all sizes 100 times 
mussels.FC.sample = NULL
for (i in 1:100) {
  s = sample(mussels.FC$Thk.length, 55, replace=TRUE)
  m = mean(s)
  mussels.FC.sample<- c(mussels.FC.sample, m)
  }
hist(mussels.FC.sample)
abline(v = 0.01994755, col= 2)

range(mussels.FC.sample) 

#Sample just 50 SH mussels that are <120 and greater than 40mm all sizes 100 times 
mussels.FC.less.120.great.40 <- mussels.FC[which(mussels.FC$Shell.Length < 120 
                                                 & mussels.FC$Shell.Length > 40),]
mussels.FC.less.120.great.40.sample = NULL
for (i in 1:100) {
  s = sample(mussels.FC.less.120.great.40$Thk.length, 55, replace=TRUE)
  m = mean(s)
  mussels.FC.less.120.great.40.sample<- c(mussels.FC.less.120.great.40.sample, m)
  }
hist(mussels.FC.less.120.great.40.sample)
abline(v = 0.01994755, col= 2)

range(mussels.FC.less.120.great.40.sample) 

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))


# Order Metadata
metadata$Site.Code <- factor(metadata$Site.Code, 
levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD","KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Remove MP from 2019 data
mussels <- mussels %>% filter(Site.Code != "MP")
# Rename BH as BMR
mussels$Site.Code[which(mussels$Site.Code == "BH")] <- "BMR"


# 2019 mussel analyses
# Summarize data
mussels_sum <- mussels %>% group_by(Site.Code) %>% summarise(mean_avg_thick = mean(Avg.Thickness))
# Bind
mussels_sum_meta <- left_join(mussels_sum, metadata_sub, by="Site.Code")

# Graph 2019 data
pdf("output/figures/Mcali_ave_thick_2019.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = mussels_sum_meta, aes(x = Longitude, y = Latitude, fill = mean_avg_thick), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 12) + ggtitle("Shell Thickness Projections")
dev.off()


# 2024 mussel analyses
# Summarize data
mussels_2024_sum <- mussels_2024 %>% group_by(Site.Code) %>% summarise(mean_STI = mean(STI), mean_integrated_thk = mean(Integrated.Thk))
# Bind
mussels_2024_sum_meta <- left_join(mussels_2024_sum, metadata, by="Site.Code")

# Graph 2024 data
pdf("output/figures/Mcali_STI_2024.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = mussels_2024_sum_meta, aes(x = Longitude, y = Latitude, fill = mean_STI), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 12) + ggtitle("STI Shell Thickness Projections")
dev.off()

# Graph 2024 data
pdf("output/figures/Mcali_integrated_thick_2024.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = mussels_2024_sum_meta, aes(x = Longitude, y = Latitude, fill = mean_integrated_thk), shape = 21, size = 5) + 
  scale_fill_gradient(low = "cyan1", high = "gray27") + 
             coord_fixed(1.3) +
  xlim(c(-128, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 12) + ggtitle("Integrated Shell Thickness Projections")
dev.off()