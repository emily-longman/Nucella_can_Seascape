# Graph MARINe data for N. canaliculata collection sites
# Note: prior to running the R script, need to load R and GDAL module on the VACC (module load R/4.4.1)

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'gameofthrones', 'vegan', 'cowplot', 'psych'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(gameofthrones)
library(vegan)
library(cowplot)
library(psych)

# ================================================================================== #

# Read in MARINe Biodiversity data
point_contact <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_point_contact_summary.csv", header=T)
quadrat <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_quadrat_summary.csv", header=T)
swath <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_swath_summary.csv", header=T)

# ================================================================================== #

# MARINe site names 
MARINe_site <- c("Fogarty Creek", "Seal Rock", "Bob Creek", "Cape Arago", "Coquille Point", "Point Saint George", 
"Shelter Cove", "Kibesillah Hill", "Windermere Point", "Bodega", "Pigeon Point", "Point Lobos", "Garrapata", 
"Point Sierra Nevada", "Piedras Blancas", "Hazards", "Stairs")

# ================================================================================== #

# Filter datasets to only include sites of interest
point_contact_filt <- point_contact %>% filter(point_contact$marine_site_name %in% MARINe_site)
quadrat_filt <- quadrat %>% filter(quadrat$marine_site_name %in% MARINe_site)
swath_filt <- swath %>% filter(swath$marine_site_name %in% MARINe_site)

# Change marine site name to factor
point_contact_filt$marine_site_name <- factor(point_contact_filt$marine_site_name, 
levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
"Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
"Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))
quadrat_filt$marine_site_name <- factor(quadrat_filt$marine_site_name, 
levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
"Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
"Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))
swath_filt$marine_site_name <- factor(swath_filt$marine_site_name, 
levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
"Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
"Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))

# ================================================================================== #

# Prey percent cover

# Summarize for multiple years of surveys
point_contact_filt_sum <- point_contact_filt %>% 
group_by(marine_site_name, latitude, longitude, georegion, state_province, species_lump) %>%
summarise(num.years=n(), mean_perc_cov = mean(percent_cover), sd_perc_cov = sd(percent_cover), se_perc_cov = sd_perc_cov/sqrt(num.years))

# Add column to specify if barnacle or mussel
point_contact_filt_sum$prey <- ifelse(point_contact_filt_sum$species_lump == "Mytilus californianus" | point_contact_filt_sum$species_lump == "Mytilus trossulus/galloprovincialis/edulis", "Mussel", "Barnacle" )

# Set colors
mycolors_prey <- c("#01a24fd5", "#487402d8", "#0d1cbc", "#006aff", "darkgreen")

# Graph prey percent cover
pdf("output/figures/GEA/enviro/MARINe/MARINe_prey_perc_cov_mycol.pdf", width = 16, height = 10)
ggplot(data = point_contact_filt_sum, aes(x=mean_perc_cov, y=marine_site_name, fill = species_lump, colour=species_lump)) + geom_point(size=2.5) + 
geom_errorbar(aes(xmin=mean_perc_cov-se_perc_cov, xmax=mean_perc_cov+se_perc_cov), width=.2) +
scale_fill_manual(values=mycolors_prey) + 
scale_colour_manual(values=mycolors_prey) +
xlab("Percent cover") + ylab("") +
theme_bw(base_size = 16) + facet_grid(. ~ prey) + 
theme(strip.background =element_rect(fill="white"))+
theme(strip.text = element_text(colour = 'black'))
dev.off()

# Graph prey percent cover (viridis colors)
pdf("output/figures/GEA/enviro/MARINe/MARINe_prey_perc_cov_viridis.pdf", width = 16, height = 10)
ggplot(data = point_contact_filt_sum, aes(x=mean_perc_cov, y=marine_site_name, fill = species_lump, colour=species_lump)) + geom_point(size=2.5) + 
geom_errorbar(aes(xmin=mean_perc_cov-se_perc_cov, xmax=mean_perc_cov+se_perc_cov), width=.2) +
scale_color_viridis(discrete=TRUE) +
xlab("Percent cover") + ylab("") +
theme_bw(base_size = 16) + facet_grid(. ~ prey) + 
theme(strip.background=element_rect(fill="white"))+
theme(strip.text = element_text(colour = 'black'))
dev.off()

# Write summary table
write.csv(point_contact_filt_sum, "data/processed/GEA/enviro_data/MARINe/point_contact_filt_sum.csv")

# ================================================================================== #

# Competitor density 

# Summarize for multiple years of surveys
quadrat_filt_sum <- quadrat_filt %>% 
group_by(marine_site_name, latitude, longitude, georegion, state_province, species_lump) %>%
summarise(num.years=n(), mean_density = mean(density_per_m2), sd_density = sd(density_per_m2), se_perc_cov = sd_density/sqrt(num.years))

mycolors_comp <- c("#5e0000", "#cf8504e3", "#ca3802dc")

# Graph competitor density 
pdf("output/figures/GEA/enviro/MARINe/MARINe_competitor_density.pdf", width = 10, height = 10)
ggplot(data = quadrat_filt_sum, aes(x=mean_density, y=marine_site_name, fill = species_lump, colour=species_lump)) + geom_point(size=3) + 
geom_errorbar(aes(xmin=mean_density-se_perc_cov, xmax=mean_density+se_perc_cov), width=.2) +
scale_fill_manual(values=mycolors_comp) + 
scale_colour_manual(values=mycolors_comp) +
xlab(bquote("Density/ m"^2)) + ylab("") +
theme_bw(base_size = 16) 
dev.off()

# Write summary table
write.csv(quadrat_filt_sum, "data/processed/GEA/enviro_data/MARINe/quadrat_filt_sum.csv")

# ================================================================================== #

# Predator density 

# Summarize for multiple years of surveys
swath_filt_sum <- swath_filt %>% 
group_by(marine_site_name, latitude, longitude, georegion, state_province) %>%
summarise(num.years=n(), mean_density = mean(density_per_m2), sd_density = sd(density_per_m2), se_perc_cov = sd_density/sqrt(num.years))

# Graph predator density 
pdf("output/figures/GEA/enviro/MARINe/MARINe_pisaster_density.pdf", width = 8, height = 10)
ggplot(data = swath_filt_sum, aes(x=mean_density, y=marine_site_name,)) + geom_point(size=3, col="#730b43") + 
geom_errorbar(aes(xmin=mean_density-se_perc_cov, xmax=mean_density+se_perc_cov), width=.2, col="#730b43") +
xlab(bquote("Pisaster Density/ m"^2)) + ylab("") +
theme_bw(base_size = 16) 
dev.off()

# Write summary table
write.csv(swath_filt_sum, "data/processed/GEA/enviro_data/MARINe/swath_filt_sum.csv")

# ================================================================================== #
# ================================================================================== #

# Use multivariate approach to look at patterns

# Load summary data
biodiversity_means <- read.csv("data/processed/GEA/enviro_data/MARINe/Biodiversity_means.csv", header=T)

biodiversity_means$marine_site_name <- factor(biodiversity_means$marine_site_name, levels=c("Fogarty Creek", "Seal Rock", "Bob Creek", "Cape Arago", "Coquille Point", "Point Saint George", 
"Shelter Cove", "Kibesillah Hill", "Windermere Point", "Bodega", "Pigeon Point", "Point Lobos", "Garrapata", 
"Point Sierra Nevada", "Piedras Blancas", "Hazards", "Stairs"))
# ================================================================================== #

# Run NMDS
set.seed(2)
biodiversity_mds <- metaMDS(biodiversity_means[,7:15], distance = "bray", trymax = 50)
biodiversity_mds

# Extract the axes of nmds and add columns with site info
data_scores <- as.data.frame(scores(biodiversity_mds, "sites"))  
data_scores$marine_site_name <- biodiversity_means$marine_site_name
data_scores$latitude <- biodiversity_means$latitude
data_scores$longitude <- biodiversity_means$longitude

#Extract the species scores
species_scores <- as.data.frame(scores(biodiversity_mds, "species"))
species_scores$species <- c("Balanus", "Chthamalus", "Mytilus_californianus", "Mytilus_spp", "Pollicipes", "N_canaliculata", "N_emar_ostrina", "N_lamellosa", "Pisaster")

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph NMDS
pdf("output/figures/GEA/enviro/MARINe/MARINe_NMDS.pdf", width = 10, height = 10)
ggplot() + 
geom_point(data=data_scores, aes(x=NMDS1, y=NMDS2, colour=marine_site_name), size=3) + 
geom_text(data=species_scores, aes(x=NMDS1, y=NMDS2, label=species), size=3) + 
coord_equal() +
theme_bw() + scale_color_manual(values=mycolors)
dev.off()

# ================================================================================== #

# Biplot of biotic data

# Perform the PCA
pca_biodiversity <- prcomp(biodiversity_means[,7:15], scale.=TRUE)

# Graph biplot
pdf("output/figures/GEA/enviro/MARINe/MARINe_biplot.pdf", width = 10, height = 10)
biplot(pca_biodiversity)
dev.off()

# Graph biplot with ggplot and ggfortify
pdf("output/figures/GEA/enviro/MARINe/MARINe_biplot_ggplot.pdf", width = 10.5, height = 8)
autoplot(pca_biodiversity, data=biodiversity_means, color="black", fill="marine_site_name", size=6, shape=21,
loadings=TRUE, loadings.label=TRUE, loadings.label.size=6) + scale_fill_manual(values=mycolors) + ylim(-0.6,0.6) + xlim(-0.6,0.6)+
theme_bw(base_size=20)
dev.off()

# ================================================================================== #

# Assess correlations among the biotic data

# Bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal
pdf("output/figures/GEA/enviro/MARINe/Correlations.pdf", width = 10, height = 10)
pairs.panels(biodiversity_means[,7:15], scale=T)
dev.off()

# Lots of them are highly correlated - need to subset, but which to chose?

###### 

# Just prey
pdf("output/figures/GEA/enviro/MARINe/Correlations_prey.pdf", width = 10, height = 10)
pairs.panels(biodiversity_means[,7:11], scale=T)
dev.off()

# Subset of prey (remove M. tross/gallo)
pdf("output/figures/GEA/enviro/MARINe/Correlations_prey_sub.pdf", width = 10, height = 10)
pairs.panels(biodiversity_means[,c(7:9,11)], scale=T)
dev.off()

######

# Just competitors and predators
pdf("output/figures/GEA/enviro/MARINe/Correlations_comp_pred.pdf", width = 10, height = 10)
pairs.panels(biodiversity_means[,12:15], scale=T)
dev.off()

# N. canaliculata and N. ostrina/emarginata are correlated - which should I drop?