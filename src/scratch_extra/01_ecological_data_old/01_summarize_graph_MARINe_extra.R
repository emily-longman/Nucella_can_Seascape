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

# Generate output directories
out_dir <- paste("data/processed/GEA/enviro_data/MARINe")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Figure directory
out_dir_fig <- paste("output/figures/enviro_data/MARINe")
if (!dir.exists(out_dir_fig)) {dir.create(out_dir_fig)}

# ================================================================================== #

# Color palettes for graphing
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
mycolors_prey <- c("#01a24fd5", "#487402d8", "#0d1cbc", "#006aff", "darkgreen")
mycolors_comp <- c("#5e0000", "#cf8504e3", "#ca3802dc")

# ================================================================================== #

# Read in MARINe Biodiversity data
point_contact <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_point_contact_summary.csv", header=T)
quadrat <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_quadrat_summary.csv", header=T)
swath <- read.csv("data/raw/MARINe/Biodiversity.Data/MARINe_Biodiversity_data_swath_summary.csv", header=T)

# ================================================================================== #

# MARINe site names that correspond to my field sites
MARINe_site <- c("Fogarty Creek", "Seal Rock", "Bob Creek", "Cape Arago", "Coquille Point", "Point Saint George", 
"Shelter Cove", "Kibesillah Hill", "Windermere Point", "Bodega", "Pigeon Point", "Point Lobos", "Garrapata", 
"Point Sierra Nevada", "Piedras Blancas", "Hazards", "Stairs")

# ================================================================================== #

# Filter datasets to only include sites of interest
point_contact_filt <- point_contact %>% filter(point_contact$marine_site_name %in% MARINe_site)
quadrat_filt <- quadrat %>% filter(quadrat$marine_site_name %in% MARINe_site)
swath_filt <- swath %>% filter(swath$marine_site_name %in% MARINe_site)

# Change MARINe site names to factors
point_contact_filt$marine_site_name <- factor(point_contact_filt$marine_site_name, 
    levels = rev(c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek")))
quadrat_filt$marine_site_name <- factor(quadrat_filt$marine_site_name, 
    levels = rev(c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek")))
swath_filt$marine_site_name <- factor(swath_filt$marine_site_name, 
    levels = rev(c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek")))

# ================================================================================== #
# ================================================================================== #

# Prey percent cover

# Remove Pollicipes polymerus
point_contact_filt <- point_contact_filt[which(point_contact_filt$species_lump !="Pollicipes polymerus"),]

# Summarize for multiple years of surveys
point_contact_filt_sum <- point_contact_filt %>% 
    group_by(marine_site_name, latitude, longitude, georegion, state_province, species_lump) %>%
    reframe(
        num_years = n(), mean = mean(percent_cover), 
        sd = sd(percent_cover), se = sd/sqrt(num_years), 
        min = range(percent_cover)[1], max = range(percent_cover)[2], harm_mean = harmonic.mean(percent_cover), 
        metric = "percent_cover") %>% as.data.frame() %>% distinct()

# Add column to specify if barnacle or mussel
point_contact_filt_sum$prey_species <- ifelse(
    point_contact_filt_sum$species_lump == "Mytilus californianus" | 
    point_contact_filt_sum$species_lump == "Mytilus trossulus/galloprovincialis/edulis", "Mussel", "Barnacle" )

# Write summary table
write.csv(point_contact_filt_sum, "data/processed/GEA/enviro_data/MARINe/point_contact_filt_sum.csv", row.names=F)

# ================================================================================== #

# Graph prey percent cover

# Temporal trend of prey
pdf("output/figures/GEA/enviro/MARINe/MARINe_prey_perc_cov_mussel_temporal_trend.pdf", width = 16, height = 10)
ggplot(data = point_contact_filt, aes(x=year, y=percent_cover, group=marine_site_name, color=marine_site_name)) + 
geom_point(size=2.5) + 
geom_line() + scale_color_manual(values=mycolors) + facet_wrap(~species_lump) +
xlab("Year") + ylab("Percent cover") +
theme_bw(base_size = 16) 
dev.off()

# Opposite order for levels
point_contact_filt_sum$marine_site_name <- factor(point_contact_filt_sum$marine_site_name, 
    levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))

# Graph summary of prey percent cover
pdf("output/figures/GEA/enviro/MARINe/MARINe_prey_perc_cov_mycol.pdf", width = 16, height = 10)
ggplot(data = point_contact_filt_sum, aes(x=mean, y=marine_site_name, fill = species_lump, colour=species_lump)) + geom_point(size=2.5) + 
geom_errorbar(aes(xmin=mean-se, xmax=mean+se), width=.2) +
scale_fill_manual(values=mycolors_prey) + 
scale_colour_manual(values=mycolors_prey) +
xlab("Percent cover") + ylab("") +
theme_bw(base_size = 16) + facet_grid(. ~ prey_species) + 
theme(strip.background =element_rect(fill="white"))+
theme(strip.text = element_text(colour = 'black'))
dev.off()

####

# Graph projections

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Cut ARA since outlier for collection
point_contact_filt_sum.17 <- point_contact_filt_sum[-which(point_contact_filt_sum$marine_site_name == "Coquille Point"),]

# Graph M. trossulus mean
pdf("output/figures/enviro/MARINe/MARINe_Mtross_17sites_map.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = point_contact_filt_sum.17[which(point_contact_filt_sum.17$species_lump == "Mytilus trossulus/galloprovincialis/edulis"),], aes(x = longitude, y = latitude, fill = mean), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  #scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name=expression(paste("Mean ", italic("M. trossulus")))) +
  scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name=NULL) +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 26) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.98, 0.52))
dev.off()
pdf("output/figures/enviro/MARINe/MARINe_Mtross_mean_17sites_map_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = point_contact_filt_sum.17[which(point_contact_filt_sum.17$species_lump == "Mytilus trossulus/galloprovincialis/edulis"),], aes(x = longitude, y = latitude, fill = mean), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  scale_fill_gradientn(colours=brewer.pal(6, "YlOrRd"), name=NULL) +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_bw(base_size = 26) + #ggtitle("Integrated Shell Thickness Projections") + 
  theme(legend.title = element_text(size = 20), legend.text = element_text(size = 20), legend.position = c(0.98, 0.52), legend.background = element_rect(color = "black", fill = "white", size = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Competitor density

# Remove Nucella lamellosa
quadrat_filt <- quadrat_filt[which(quadrat_filt$species_lump !="Nucella lamellosa"),]

# Summarize for multiple years of surveys
quadrat_filt_sum <- quadrat_filt %>% 
    group_by(marine_site_name, latitude, longitude, georegion, state_province, species_lump) %>%
    reframe(num_years=n(), mean = mean(density_per_m2), 
    sd = sd(density_per_m2), se = sd/sqrt(num_years), 
    min = range(density_per_m2)[1], max = range(density_per_m2)[2], harm_mean = harmonic.mean(density_per_m2), 
    metric = "density") %>% as.data.frame() %>% distinct()

# Write summary table
write.csv(quadrat_filt_sum, "data/processed/GEA/enviro_data/MARINe/quadrat_filt_sum.csv", row.names=F)

# ================================================================================== #

# Graph competitor density

# Opposite order for levels
quadrat_filt_sum$marine_site_name <- factor(quadrat_filt_sum$marine_site_name, 
    levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))

# Graph summary of competitor density
pdf("output/figures/GEA/enviro/MARINe/MARINe_competitor_density.pdf", width = 10, height = 10)
ggplot(data = quadrat_filt_sum, aes(x=mean, y=marine_site_name, fill = species_lump, colour=species_lump)) + geom_point(size=3) + 
geom_errorbar(aes(xmin=mean-se, xmax=mean+se), width=.2) +
scale_fill_manual(values=mycolors_comp) +
scale_colour_manual(values=mycolors_comp) +
xlab(bquote("Density/ m"^2)) + ylab("") +
theme_bw(base_size = 16)
dev.off()

# Temporal trend of competitors
pdf("output/figures/GEA/enviro/MARINe/MARINe_competitors_density_temporal_trend.pdf", width = 16, height = 8)
ggplot(data = quadrat_filt, aes(x=year, y=density_per_m2, group=marine_site_name, color=marine_site_name)) + 
geom_point(size=2.5) + 
geom_line() + scale_color_manual(values=mycolors) + facet_wrap(~species_lump) +
xlab("Year") + ylab("Density per m2") +
theme_bw(base_size = 16) 
dev.off()

# ================================================================================== #
# ================================================================================== #

# Predator (pisaster) density

# Summarize for multiple years of surveys
swath_filt_sum <- swath_filt %>% 
    group_by(marine_site_name, latitude, longitude, georegion, state_province, species_lump) %>%
    reframe(num_years=n(), mean = mean(density_per_m2), 
    sd = sd(density_per_m2), se = sd/sqrt(num_years), 
    min = range(density_per_m2)[1], max = range(density_per_m2)[2], harm_mean = harmonic.mean(density_per_m2), 
    metric = "density") %>% as.data.frame() %>% distinct()

# Write summary table
write.csv(swath_filt_sum, "data/processed/GEA/enviro_data/MARINe/swath_filt_sum.csv", row.names=F)

# ================================================================================== #

# Graph predator (Pisaster) density

# Temporal trend of Pisaster
pdf("output/figures/GEA/enviro/MARINe/MARINe_pisaster_density_temporal_trend.pdf", width = 16, height = 10)
ggplot(data = swath_filt, aes(x=year, y=density_per_m2, group=marine_site_name, color=marine_site_name)) + geom_point(size=2.5) + 
geom_line() + scale_color_manual(values=mycolors) +
xlab("Year") + ylab("Pisaster density") +
theme_bw(base_size = 16)
dev.off()

# Opposite order for levels
swath_filt_sum$marine_site_name <- factor(swath_filt_sum$marine_site_name, 
    levels = c("Stairs", "Hazards", "Piedras Blancas", "Point Sierra Nevada", "Garrapata", "Point Lobos", 
    "Pigeon Point","Bodega","Windermere Point", "Kibesillah Hill", "Shelter Cove", "Point Saint George", 
    "Coquille Point", "Cape Arago", "Bob Creek", "Seal Rock", "Fogarty Creek"))

# Graph Pisaster density
pdf("output/figures/GEA/enviro/MARINe/MARINe_pisaster_density.pdf", width = 8, height = 10)
ggplot(data = swath_filt_sum, aes(x=mean, y=marine_site_name,)) + geom_point(size=3, col="#730b43") + 
geom_errorbar(aes(xmin=mean-se, xmax=mean+se), width=.2, col="#730b43") +
xlab(bquote("Pisaster Density/ m"^2)) + ylab("") +
theme_bw(base_size = 16)
dev.off()

# ================================================================================== #
# ================================================================================== #

# Join the datasets

# Rename species lump for quadrat_filt_cong_2_sum
colnames(quadrat_filt_sum)[6] <- "species_lump"

# Join the 3 datasets - remove prey_species column from point_contact_filt_sum
marine_sum <- rbind(point_contact_filt_sum[,-15], quadrat_filt_sum, swath_filt_sum)

# Write summary table
write.csv(marine_sum, "data/processed/GEA/enviro_data/MARINe/marine_long_sum.csv", row.names=F)

# ================================================================================== #
# ================================================================================== #

# Use multivariate approach to look at patterns

# Note - will need to reformat the data if plan to use these

# Load summary data
#biodiversity_means <- read.csv("data/processed/GEA/enviro_data/MARINe/Biodiversity_means.csv", header=T)

#biodiversity_means$marine_site_name <- factor(biodiversity_means$marine_site_name, levels=c("Fogarty Creek", "Seal Rock", "Bob Creek", "Cape Arago", "Coquille Point", "Point Saint George", 
#"Shelter Cove", "Kibesillah Hill", "Windermere Point", "Bodega", "Pigeon Point", "Point Lobos", "Garrapata", 
#"Point Sierra Nevada", "Piedras Blancas", "Hazards", "Stairs"))

# ================================================================================== #

# Run NMDS
set.seed(2)
biodiversity_mds <- metaMDS(marine_sum_means[,6:14], distance = "bray", trymax = 50)
biodiversity_mds

# Extract the axes of nmds and add columns with site info
data_scores <- as.data.frame(scores(biodiversity_mds, "sites"))
data_scores$marine_site_name <- marine_sum_means$marine_site_name
data_scores$latitude <- marine_sum_means$latitude
data_scores$longitude <- marine_sum_means$longitude

#Extract the species scores
species_scores <- as.data.frame(scores(biodiversity_mds, "species"))
species_scores$species <- c("Balanus", "Chthamalus", "Mytilus_californianus", "Mytilus_spp", "Pollicipes", "N_canaliculata", "N_emar_ostrina", "N_lamellosa", "Pisaster")

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
pca_biodiversity <- prcomp(marine_sum_means[,6:14], scale.=TRUE)

# Graph biplot
pdf("output/figures/GEA/enviro/MARINe/MARINe_biplot.pdf", width = 10, height = 10)
biplot(pca_biodiversity)
dev.off()

# Graph biplot with ggplot and ggfortify
pdf("output/figures/GEA/enviro/MARINe/MARINe_biplot_ggplot.pdf", width = 10.5, height = 8)
autoplot(pca_biodiversity, data=marine_sum_means, color="black", fill="marine_site_name", size=6, shape=21,
loadings=TRUE, loadings.label=TRUE, loadings.label.size=6) + scale_fill_manual(values=mycolors) + ylim(-0.6,0.6) + xlim(-0.6,0.6)+
theme_bw(base_size=20)
dev.off()

# ================================================================================== #

# Assess correlations among the biotic data

# Bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal
pdf("output/figures/GEA/enviro/MARINe/Correlations.pdf", width = 10, height = 10)
pairs.panels(marine_sum_means[,6:14], scale=T)
dev.off()

# Lots of them are highly correlated - need to subset, but which to chose?