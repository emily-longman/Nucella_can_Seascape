# Compare correlation coef between SNPs in gene g27343 and those with similar AF that are not on that scaffold

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis', 'colorspace'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)
library(colorspace)

# ================================================================================== #

# Set seed
set.seed(1234)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/outliers")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load Baypass alpha(ij)
load("data/processed/baypass/abiotic/baypass.ph.sum.RData")

# Load AFs for SNPs that beat POD in gene g27343
load("data/processed/baypass/afs.ph.g27343.BF.POD.RData")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)
# Extract just mean pH and rename location to Site
ph <- bio_oracle_sites_2010[,c(1,2,3,11,13)]
ph <- ph %>% rename(Site = location)

# Load PCA data for demography
pca.df <- read.csv("data/processed/outlier_analyses/pca.csv")
colnames(pca.df)[1] <- "Site"
PC1 <- pca.df[,c(1,2)]

# ================================================================================== #

# Extract data for just outlier SNPs on g27343

# Join outlier AF dataset with baypass corrected AF
ph.g27343 <- left_join(afs.ph.g27343.BF.POD, baypass.ph.sum, by = c("SNP_id", "Site"))
# Join with pH data
ph.g27343 <- left_join(ph.g27343, ph)

# ================================================================================== #

# Summarize across sites to calculate global AFs
baypass.ph.sum.ag <- baypass.ph.sum %>% group_by(MRK, chr, pos, allele1, allele2, SNP_id) %>%
    reframe(mean_M_P=mean(M_P.mean, na.rm=T))

#####

# Identify 1,000 SNPs with similar mean AF to those outlier SNPs in g27343

# Extract summary AF data for all outlier SNPs in g27343
baypass.ph.sum.ag.g27343 <- baypass.ph.sum.ag %>% filter(SNP_id %in% unique(afs.ph.g27343.BF.POD$SNP_id))
# Calc summary stats
baypass.g27343.mean.MP <- mean(baypass.ph.sum.ag.g27343$mean_M_P)
baypass.g27343.sd.MP <- sd(baypass.ph.sum.ag.g27343$mean_M_P)
baypass.g27343.range.MP <- range(baypass.ph.sum.ag.g27343$mean_M_P)

# Extract a sample of 1000 SNPs that aren't on that contig and have a mean AF within one sd of the mean AF for g27343
baypass.ph.sum.ag.sample <- baypass.ph.sum.ag %>% 
    filter(SNP_id != unique(baypass.ph.sum.ag.g27343$chr)) %>% 
    filter(mean_M_P > baypass.g27343.mean.MP-baypass.g27343.sd.MP & mean_M_P < baypass.g27343.mean.MP+baypass.g27343.sd.MP) |> 
    slice_sample(n = 1000)

# ================================================================================== #

# Extract per site AF for each of the 1,000 randomly sampled SNPs
ph.sample <- baypass.ph.sum %>% filter(SNP_id %in% baypass.ph.sum.ag.sample$SNP_id)

# Join with pH data
ph.sample <- left_join(ph.sample, ph, by="Site")

# ================================================================================== #

# Join with PC1 data
ph.g27343.PC1 <- left_join(ph.g27343, PC1)
# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
baypass.cor.g27343.PC1.lat <- foreach(i=1:length(unique(ph.g27343.PC1$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = ph.g27343.PC1 %>% filter(SNP_id == unique(ph.g27343.PC1$SNP_id)[i])
    # Model
    mod.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + ph_mean, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(ph.g27343.PC1$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "outliers")
}
# Create column with absolute value of cor
baypass.cor.g27343.PC1.lat$abs.estimate <- abs(baypass.cor.g27343.PC1.lat$estimate)

# Join with PC1 data
sample.PC1 <- left_join(ph.sample, PC1)
# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
baypass.cor.sample.PC1.lat <- foreach(i=1:length(unique(sample.PC1$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = sample.PC1 %>% filter(SNP_id == unique(sample.PC1$SNP_id)[i])
    # Model
    mod.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + ph_mean, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(sample.PC1$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "sample")
}
# Create column with absolute value of cor
baypass.cor.sample.PC1.lat$abs.estimate <- abs(baypass.cor.sample.PC1.lat$estimate)

# ================================================================================== #

# Join - PC1 and lat incorporated
baypass.cor.g27343.PC1.lat.join <- rbind(baypass.cor.g27343.PC1.lat, baypass.cor.sample.PC1.lat)

# Colors
cols <- c(head(hcl.colors(12, palette = "BrwnYl"), 1), "#757474")
cols_alt <- c("#9c52596d", "#757474")

# Graph
pdf("output/figures/baypass/outliers/Correlation_abs_ph_baypass_PC1_lat_density_fill.pdf", width = 11, height = 4)
ggplot(baypass.cor.g27343.PC1.lat.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,0.5) +
  scale_fill_manual(values = cols) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30)
dev.off()

# Graph alt
pdf("output/figures/baypass/outliers/Correlation_abs_ph_baypass_PC1_lat_density_fill_alt.pdf", width = 6, height = 6)
ggplot(baypass.cor.g27343.PC1.lat.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,0.478) +
  scale_fill_manual(values = cols_alt) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

