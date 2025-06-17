# Format Bio-Oracle data (https://www.bio-oracle.org/index.php)

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'ggfortify'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggfortify)

# ================================================================================== #


# Read in Bio-oracle 2010 data for sites
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Set location as leveled factor
bio_oracle_sites_2010$location <- factor(bio_oracle_sites_2010$location, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# ================================================================================== #

# Biplot of environmental data

# Perform the PCA
pca_bio_oracle_sites_2010 <- prcomp(bio_oracle_sites_2010[,4:12], scale.=TRUE)

# Set color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph biplot
pdf("output/figures/GEA/enviro/Bio-oracle/Bio-oracle_biplot.pdf", width = 10, height = 10)
biplot(pca_bio_oracle_sites_2010)
dev.off()

# Graph biplot with ggplot and ggfortify
pdf("output/figures/GEA/enviro/Bio-oracle/Bio-oracle_biplot.pdf", width = 11, height = 10)
autoplot(pca_bio_oracle_sites_2010, data=bio_oracle_sites_2010, color="black", fill="location", size=6, shape=21,
loadings=TRUE, loadings.label=TRUE, loadings.label.size=6) + scale_fill_manual(values=mycolors) + ylim(-0.48,0.48) + xlim(-0.48,0.48)+
theme_bw(base_size=20)
dev.off()

