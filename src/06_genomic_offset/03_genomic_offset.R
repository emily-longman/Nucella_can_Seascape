# Perform genomic offset analyses with Bio-Oracle data
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

# Clear memory
rm(list=ls()) 
# Global setting to stop scientific notation
options(scipen = 999)

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'terra', 'raster', 'SeqArray', 'poolfstat', 'vegan', 'cowplot', 'psych'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(terra)
library(raster)
require(poolfstat)
require(vegan)
library(cowplot)
library(psych)

# Load gradientForest
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
require(gradientForest)

# ================================================================================== #

# Generate output directories

# Genomic offset data directory
out_dir <- paste("data/processed/genomic_offset")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# Genomic offset figure directory
out_fig_dir <- paste("output/figures/genomic_offset")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load SNPs of interest (baypass POD outlier SNPs - 8914 SNPs SNPs)
load("data/processed/genomic_offset/NC.xtx.pos.Rdata")
# Create SNP_id column for outlier SNP list
NC.xtx.pos <- NC.xtx.pos %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# Bio-oracle data
# Load bio-oracle environmental data - present (2000-2010)
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)
# Load bio-oracle environmental data - future (2080-2090)
bio_oracle_sites_ssp585_2090 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_ssp585_2090.csv")

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# ================================================================================== #

# For environmental data need a Site-by-Enviro matrix

# Rename "location" column as "sampleId"
names(bio_oracle_sites_2010)[names(bio_oracle_sites_2010) == "location"] <- "Site"

# Rename rows
row.names(bio_oracle_sites_2010) <- bio_oracle_sites_2010$Site
# Subset bio-oracle data so just environmental variables
bio_oracle_sites_2010_sub <- bio_oracle_sites_2010[, c(4:12)]

# Rename rows of future data as sites
rownames(bio_oracle_sites_ssp585_2090) <- rownames(bio_oracle_sites_2010_sub)
# Subset bio-oracle data so just environmental variables
bio_oracle_sites_ssp585_2090_sub <- bio_oracle_sites_ssp585_2090[, c(4:12)]

# ================================================================================== #

# Select only non-correlated variables - based on biplot

# Bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal
pdf("output/figures/genomic_offset/Correlations.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_sub, scale=T)
dev.off()

# Keep variables that are most relevant based on model enrichment
bio_oracle_sites_2010_notcorr <- bio_oracle_sites_2010_sub[,c(1,5,8)]

# ================================================================================== #
# ================================================================================== #

# Calculate allele frequencies

# Subset pooldata for Baypass outlier SNPs

# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  snp.info

# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Filter pooldata for Baypass SNPs of interest
selected_SNPs <- snp.info %>% filter(snp.info$SNP_id %in% NC.xtx.pos$SNP_id)
# Get index of SNPs
selected_SNPs_index <- as.integer(selected_SNPs$rs.id)

# Subset the pooldata object using the selected SNP indices
pooldata_subset <- pooldata.subset(pooldata, snp.index = selected_SNPs_index)

# ================================================================================== #

# Extract and manipulate snp info for significant SNPs

# Extract SNP info for significant SNPs
pooldata_subset@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  snp.info.subset

# Rename columns
names(snp.info.subset)[1:2] = c("chr","pos")
# Make snp_id column
snp.info.subset %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info.subset

# ================================================================================== #

# Extract and manipulate coverage for significant SNPs

# Extract read count and coverage data for SNPs
ref_count <- pooldata_subset@refallele.readcount
coverage <- pooldata_subset@readcoverage

# Change to data frame
coverage %>% as.data.frame -> cov
# Rename columns (19 sites)
names(cov) = c(pooldata_subset@poolnames)

# Transpose coverage
cov_t <- t(cov)

# ================================================================================== #

# Calculate allele frequencies for significant SNPs

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage

# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata_subset@poolnames)

# Transpose allele freqs
afs_t <- t(afs)

# Rename columns as SNP id
colnames(afs_t) <- snp.info.subset$SNP_id

# ================================================================================== #

# Sample size of each pool
nSnail=20

# Calculate mean effective coverage
nEff <- round((cov_t*2*nSnail)/(cov_t+2*nSnail-1))
# Calculate the effective allele freq
af_nEff <- round((afs_t*nEff)/nEff)

# Check for NAs -- gradient forest can't run if the data includes NAs
which(is.na(af_nEff))

# ================================================================================== #
# ================================================================================== #

# Run gradient forest

# Check size of input data
nSites <- dim(bio_oracle_sites_2010_sub)[1] #19
nSpec <- dim(af_nEff)[2] #8914

# Calc maximum number of splits for conditional permutation
maxLevel <- log2(nSites * 0.368/2)
# Note: For conditional permutation, the predictor to be assessed is permuted only within blocks of the dataset defined by splits in the given tree
# on any other predictors correlated above a certain threshold and up to a maximum number of splits set by the maxLevel option 

# Run gradient forest
gf <- gradientForest(cbind(af_nEff, bio_oracle_sites_2010_sub), predictor.vars = colnames(bio_oracle_sites_2010_sub),
                                  response.vars = colnames(af_nEff), ntree = 5000, 
                                  maxLevel = maxLevel, trace = T, corr.threshold = 0.5)
# Note: Filter by correlation threshold of 0.5
## Warnings:
# In randomForest.default(x = X, y = spec_vec, maxLevel = maxLevel,  ... :
# The response has five or fewer unique values.  Are you sure you want to do regression?


gf
# Important variables:
# [1] thetao_min   thetao_range thetao_mean  thetao_max   ph_mean
#UPDATED
# [1] thetao_min   thetao_range thetao_mean  o2_mean      thetao_max  

# ================================================================================== #

# Graph Gradient Forest Results

# Graph predictor importance (This show the mean accuracy importance and the mean importance weighted by species R2)
pdf("output/figures/genomic_offset/predict_importance.pdf", width = 8, height = 8)
plot(gf, plot.type = "O")
dev.off()

# Extract most important variables
most_important <- names(importance(gf))[1:5]

# Splits density plot (This shows binned split importance and location on each gradient (spikes), kernel density of splits (black lines), of observations
# (red lines) and of splits standardised by observations density (blue lines). These show where important changes in the abundance of multiple species are occurring along the gradient
pdf("output/figures/genomic_offset/splits_density.pdf", width = 8, height = 8)
plot(gf, plot.type = "S", imp.vars = most_important, leg.posn = "topleft", cex.legend = 0.8, cex.axis = 0.6,
cex.lab = 1, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
dev.off()
# Warning messages:
# 1: In density.default(splits, weight = w/sum(w), from = rX[1], to = rX[2]) :
#  Selecting bandwidth *not* using 'weights'

# Species cumulative plot (For each species shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations)
pdf("output/figures/genomic_offset/species_cumulative_plot.pdf", width = 8, height = 8)
plot(gf, plot.type = "C", imp.vars = most_important, show.overall = F, legend = T, leg.posn = "topleft",
leg.nspecies = 5, cex.lab = 1.5, cex.legend = 0.4, cex.axis = 1, line.ylab = 0.9, 
par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
dev.off()

# Predictor cumulative plot (common.scale=T ensures that plots for all predictors have the same y-scale)
# (For each predictor shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species.)
pdf("output/figures/genomic_offset/predictor_cumulative_plot.pdf", width = 8, height = 8)
plot(gf, plot.type = "C", imp.vars = most_important, show.species = F, common.scale = T, 
cex.axis = 1, cex.lab = 1.5, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
dev.off()

# Fit of Random Forest
pdf("output/figures/genomic_offset/grad_forest_fit.pdf", width = 8, height = 8)
plot(gf, plot.type = "P", show.names = T, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)
dev.off()

# ================================================================================== #

# Extract gradient forest results and graph accuracy importance

# Accuracy importance is a measure of how much worse the model gets (in terms of mean square prediction error, MSPE) 
# when the values of a given predictor are randomly permuted among the training data

# Extract overall importance - (accuracy importance)
overall.imp <- as.data.frame(gf$overall.imp)
# Rename column
names(overall.imp)[1] <- "Importance"
# Re-format and add variables as a column rather than rownames
overall.imp <- tibble::rownames_to_column(overall.imp, "Variable")

# Re-graph accuracy importance
pdf("output/figures/genomic_offset/accuracy_importance_ggplot.pdf", width = 8, height = 8)
ggplot(data=overall.imp, aes(x=reorder(Variable, Importance),y=(Importance))) + 
geom_bar(stat="identity", colour="black", width=0.8) +
theme_bw() +
geom_hline(yintercept = 0) +
coord_flip() +
theme(panel.grid.major = element_blank(), legend.position.inside=c(0.8,0.3),legend.key.size = unit(1, 'lines'),
        legend.title = element_blank(), panel.border = element_blank(),axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12), axis.line.x = element_line(colour="black"), 
        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),axis.text.y=element_text(size = 12),axis.title=element_text(size=14,color="black"),
        panel.background = element_blank())+ylab("Accuracy Importance")+ xlab("")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Gradient Forest Predictions - using entire model

# Transform present data with gradient forest model
predOUT_present <- predict(gf, bio_oracle_sites_2010_sub)

# Transform future data with gradient forest model
predOUT_2090 <- predict(gf, bio_oracle_sites_ssp585_2090_sub)

# ================================================================================== #

# Calculate Offset

# Calculate offset for each environmental variable (predOUT_2090 - predOUT_present)
df_diff_squared <- data.frame(
  lapply(1:ncol(predOUT_2090), function(i) {
    (predOUT_2090[,i] - predOUT_present[,i])^2
  }))
# Rename columns
colnames(df_diff_squared) <- colnames(predOUT_2090)

# Extract r squared
check_rsq <- gf$res

# Summarize R2
# Mean R2 of the individual SNP models can serve as a measure of deviance explained
mean(check_rsq$rsq)*100 #13.19323
range(check_rsq$rsq)*100 #0.00535867 66.29889841

# ================================================================================== #

# Change data to long format

# Add column with site names
df_diff_squared$Site <- rownames(predOUT_2090)

# Change to long format
df_long <- df_diff_squared %>%
  pivot_longer(cols = -Site, names_to = "Variable", values_to = "Value")

# ================================================================================== #

# Graph offset by site and environmental variable

# Rename variables
df_long <- df_long %>% 
  mutate(
    Variable = case_when(
      grepl("chl_mean", Variable)   ~ "Chlorophyll (mean)",
      grepl("o2_mean", Variable)  ~ "O2 (mean)",
      grepl("ph_mean", Variable)   ~ "pH (mean)",
      grepl("ph_min", Variable)~ "pH (min)",
      grepl("so_mean", Variable)~ "Salinity (mean)",
      grepl("thetao_max", Variable)~ "Temperature (max)",
      grepl("thetao_mean", Variable)~ "Temperature (mean)",
      grepl("thetao_min", Variable)~ "Temperature (min)",
      grepl("thetao_range", Variable)~ "Temperature (range)",
    )
  )

# Change Site to a factor and order N to S 
df_long$Site <- factor(df_long$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Change Enviro var to a factor and order by importance
df_long$Variable <- factor(df_long$Variable, levels=c("Temperature (min)", "Temperature (mean)", "Temperature (range)", "Temperature (max)", 
"pH (mean)", "O2 (mean)", "Salinity (mean)", "Chlorophyll (mean)", "pH (min)"))

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph offset
pdf("output/figures/genomic_offset/genomic_offset.pdf", width = 12, height = 8)
ggplot(df_long, aes(x = Site, y = Value, fill = Site)) + geom_col() +  
facet_wrap(~ Variable, scales = "free_y") +
scale_fill_manual(values = mycolors) +
xlab("Site") + ylab("Genomic offset") +
theme_minimal(base_size=16) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
theme(strip.text = element_text(face = "bold", size = 16)) +
guides(fill = "none")
dev.off()

# Change color pallet
viridiscolors <- viridis(n=19)
# Graph offset
pdf("output/figures/genomic_offset/genomic_offset_viridis.pdf", width = 12, height = 8)
ggplot(df_long, aes(x = Site, y = Value, fill = Site)) + geom_col() +  
facet_wrap(~ Variable, scales = "free_y") +
scale_fill_manual(values = viridiscolors) +
xlab("Site") + ylab("Genomic offset") +
theme_minimal(base_size=16) +
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
theme(strip.text = element_text(face = "bold", size = 16)) +
guides(fill = "none")
dev.off()

# ================================================================================== #

# Graph Overall Genomic Offset

# Sum all columns to create overall offset values
df_diff_squared$sum_all_columns <- rowSums(df_diff_squared[1:9])

# Change Site to a factor and order N to S 
df_diff_squared$Site <- factor(df_diff_squared$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph offset
pdf("output/figures/genomic_offset/genomic_offset_diff_sq.pdf", width = 12, height = 8)
ggplot(df_diff_squared, aes(x = Site, y = sum_all_columns, fill = Site)) +
  geom_col() +  # Use geom_col() for bar plot, or geom_point() for scatter
  theme_minimal(base_size=16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Site", y = "Genomic offset")+
  scale_fill_manual(values = mycolors) + guides(fill = "none")
  dev.off()
pdf("output/figures/genomic_offset/genomic_offset_diff_sq_viridis.pdf", width = 12, height = 8)
ggplot(df_diff_squared, aes(x = Site, y = sum_all_columns, fill = Site)) +
  geom_col() +  # Use geom_col() for bar plot, or geom_point() for scatter
  theme_minimal(base_size=16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Site", y = "Genomic offset")+
  scale_fill_manual(values = viridiscolors) + guides(fill = "none")
  dev.off()

# ================================================================================== #

# Load Raster Maps

# Load tif files generated from biooracle
env_tif_present <- list.files("data/processed/GEA/enviro_data/Bio-oracle/tif_files/", pattern ="bio-oracle_present_", full.names = TRUE)
env_tif_ssp585 <- list.files("data/processed/GEA/enviro_data/Bio-oracle/tif_files/", pattern ="bio-oracle_ssp585", full.names = TRUE)

# Stack tif files
env_stack_present <- stack(env_tif_present)
env_stack_ssp585 <- stack(env_tif_ssp585)

# Rename layers in stack
names(env_stack_present) <- sub("bio.oracle_present_", "", names(env_stack_present))
names(env_stack_present) <- sub("_raster", "", names(env_stack_present))
names(env_stack_ssp585) <- sub("bio.oracle_ssp585", "", names(env_stack_ssp585))
names(env_stack_ssp585) <- sub("_raster", "", names(env_stack_ssp585))

# ================================================================================== #

# Graphing functions

# PCA to raster
pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  # Assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  # Scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  # Assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  # Stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned

RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}

# ================================================================================== #

# Create empty raster template

# Set geographic constraints
latitude_range <- c(34, 45)
longitude_range <- c(-125, -120)
# Set study extent
study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])
# Raster resolution (focal cells of bio-oracle data are at 0.05 degree resolution)
raster_resolution <- 0.05

# Create raster layer object - specify study extent, resolution and coordinate reference system
study_raster <- raster(study_extent, res=raster_resolution, crs="+proj=longlat +datum=WGS84")

# Fill dataset with NA
values(study_raster) <- NA

# Graph empty template
pdf("output/figures/genomic_offset/Raster_template.pdf", width = 5, height = 5)
plot(study_raster, main = "West Coast Raster Template")
dev.off()

# Rename study raster to "mask" to match genomic offset functions
mask <- study_raster

# ================================================================================== #

# Format data for graphing functions

# Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
# (2) a raster mask of the study region to which the RGB data will be written

# Change raster stack to data frame with lat and long
env_df <- as.data.frame(env_stack_present, xy = TRUE, cell = TRUE)
env_df_fut <- as.data.frame(env_stack_ssp585, xy = TRUE, cell = TRUE)

# Create cell ID column
env_df$cell <- cellFromXY(env_stack_present, env_df[, c("x", "y")])
env_df_fut$cell <- cellFromXY(env_stack_ssp585, env_df_fut[, c("x", "y")])

# Reorder columns to place "cell" first to match format in graphing functions
env_df <- env_df %>% dplyr::select(cell, everything())
env_df_fut <- env_df_fut %>% dplyr::select(cell, everything())

# Omit missing data
env_df <- na.omit(env_df)
env_df_fut <- na.omit(env_df_fut)

# Select just the environmental variables for graphing functions
env_df_sub <- env_df %>% dplyr::select("chl_mean":"thetao_range")
env_df_fut_sub <- env_df_fut %>% dplyr::select("chl_mean":"thetao_range")

# ================================================================================== #

# Transform env using gf models
predRefind <- predict(gf, env_df_sub)
predRefindfut <- predict(gf, env_df_fut_sub)

# ================================================================================== #

# Graph maps

# Present data map
present_RGBmap <- pcaToRaster(predRefind, mask, env_df$cell)
# Graph
pdf("output/figures/genomic_offset/present_RGBmap.pdf", width = 8, height = 6)
plotRGB(present_RGBmap)
dev.off()
# Write Raster
writeRaster(present_RGBmap, filename = "output/figures/genomic_offset/present_RGBmap.tif", format = "GTiff")

# Future data map
future_RGBmap <- pcaToRaster(predRefindfut, mask, env_df_fut$cell)
# Graph
pdf("output/figures/genomic_offset/future_RGBmap.pdf", width = 8, height = 6)
plotRGB(future_RGBmap)
dev.off()
# Write Raster
writeRaster(future_RGBmap, filename = "output/figures/genomic_offset/future_RGBmap.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (future and present)
diff_RGBmap <- RGBdiffMap(predRefind, predRefindfut, rast=mask, mapCells=env_df$cell)
# Graph
pdf("output/figures/genomic_offset/diff_raster.pdf", width = 8, height = 6)
plot(diff_RGBmap[[2]])
dev.off()
# Write Raster
writeRaster(diff_RGBmap[[2]], filename = "output/figures/genomic_offset/diff_RGBmap.tif", format="GTiff", overwrite=TRUE)

# ================================================================================== #

# Calculate and map genomic offset 

# Script assumes:
# (1) a dataframe of transformed env. variables for CURRENT climate.
# (2) a dataframe named env_trns_future containing extracted raster data of 
# env. variables for FUTURE a climate scenario, same structure as env_trns

# Calculate euclidean distance between current and future genetic spaces  
genOffset <- sqrt((predRefindfut[,1]-predRefind[,1])^2+(predRefindfut[,2]-predRefind[,2])^2
                     +(predRefindfut[,3]-predRefind[,3])^2+(predRefindfut[,4]-predRefind[,4])^2
                     +(predRefindfut[,5]-predRefind[,5])^2+(predRefindfut[,6]-predRefind[,6])^2
                     +(predRefindfut[,7]-predRefind[,7])^2)


# Assign gen offset to mask - can be tricky if current/future climate rasters are not identical in terms of # cells, extent, etc.
mask[env_df_fut$cell] <- genOffset

# Graph genomic offset
pdf("output/figures/genomic_offset/genomic_offset_map.pdf", width = 8, height = 6)
plot(mask)
dev.off()

# ================================================================================== #

# Graph map with ggplot

# Convert raster to a dataframe for ggplot
mask_df <- as.data.frame(rasterToPoints(mask))
colnames(mask_df) <- c("Longitude", "Latitude", "Genomic_Offset")

# Graph with ggplot2
pdf("output/figures/genomic_offset/genomic_offset_ggplot.pdf", width = 8, height = 6)
ggplot(mask_df, aes(x = Longitude, y = Latitude, fill = Genomic_Offset)) +
  geom_raster() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_gradientn(colours=brewer.pal(5, "RdBu")) +
  scale_fill_gradientn(colours=brewer.pal(10, "RdGy")) +
  #scale_fill_viridis(option = "inferno", na.value = "transparent", name = "Genomic Offset") + 
  theme_minimal(base_size = 12) +
  coord_fixed() +  # Keeps aspect ratio
  xlab("Longitude") + ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), axis.title = element_text(size = 16), axis.text = element_text(size = 16))
dev.off()

# Make better resolution - increase resolution to 2x
mask_highres <- disaggregate(mask, fact = 2) 
mask_highres <- as.data.frame(rasterToPoints(mask_highres))
colnames(mask_highres) <- c("Longitude", "Latitude", "Genomic_Offset")

# Graph with ggplot2
# Note: Use geom_tile instead of geom_raster for higher res
pdf("output/figures/genomic_offset/genomic_offset_ggplot2_brewer.pdf", width = 8, height = 6)
ggplot(mask_highres, aes(x = Longitude, y = Latitude, fill = Genomic_Offset)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdGy") +
  theme_minimal() +
  coord_fixed() +
  xlab("Longitude") + ylab("Latitude") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 10), axis.text = element_text(size = 10))
dev.off()



