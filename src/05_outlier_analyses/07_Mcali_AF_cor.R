# Compare correlation of Baypass alpha with M. californianus cross sectional thickness compared to random sample

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis', 'stats'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)
library(stats)

# ================================================================================== #

# Set seed
set.seed(1234)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/outliers")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# M. cali thk
load("data/processed/baypass/biotic/bf.Mcali.IntThk.sum.Rdata")

# Load Baypass alpha(ij)
load("data/processed/baypass/biotic/baypass.Mcali.sum.RData")

# Load M californianus shell thickness data
Mcalifornianus_data <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean_18pop.csv", header=T)
Mcali <- Mcalifornianus_data[,c(1,3,4,7)]
Mcali <- Mcali %>% rename(Site = Site.Code, latitude = Latitude, longitude = Longitude)

# Load PCA data for demography
pca.df <- read.csv("data/processed/outlier_analyses/pca.csv")
colnames(pca.df)[1] <- "Site"
PC1 <- pca.df[,c(1,2)]

# ================================================================================== #

# Extract outlier SNP

# Identify top SNP
Mcali.top.outlier <- Mcali.cohens.fsq.PC1.lat.bf[which(Mcali.cohens.fsq.PC1.lat.bf$cohens.fsq > 8 & Mcali.cohens.fsq.PC1.lat.bf$bf_db.mean >= 25),]

# Get baypass alpha for top SNP
alpha.Mcali.top.outlier <- baypass.Mcali.sum %>% filter(SNP_id == Mcali.top.outlier$SNP_id)

# Join with Mcali data
alpha.Mcali.top.outlier <- left_join(alpha.Mcali.top.outlier, Mcali)
# Join with PC1 data
alpha.Mcali.top.outlier <- left_join(alpha.Mcali.top.outlier, PC1)

# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
baypass.cor.Mcali.out.PC1.lat <- foreach(i=1:length(unique(alpha.Mcali.top.outlier$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = alpha.Mcali.top.outlier %>% filter(SNP_id == unique(alpha.Mcali.top.outlier$SNP_id)[i])
    # Model
    mod.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + mean_integrated_thk, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(alpha.Mcali.top.outlier$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "outliers")
}
# Create column with absolute value of cor
baypass.cor.Mcali.out.PC1.lat$abs.estimate <- abs(baypass.cor.Mcali.out.PC1.lat$estimate)

# ================================================================================== #

# Summarize across sites to calculate global AFs
baypass.Mcali.sum.ag <- baypass.Mcali.sum %>% group_by(MRK, chr, pos, allele1, allele2, SNP_id) %>%
    reframe(mean_M_P=mean(M_P.mean, na.rm=T))

#####

# Identify 1,000 SNPs with similar mean AF to those outlier SNP

# Extract summary AF data for outlier SNP
baypass.Mcali.sum.ag.out <- baypass.Mcali.sum.ag %>% filter(SNP_id %in% unique(alpha.Mcali.top.outlier$SNP_id))
baypass.Mcali.mean.MP <- baypass.Mcali.sum.ag.out$mean_M_P

# Extract a sample of 1000 SNPs that aren't on that contig and have a mean alpha within 0.15
baypass.Mcali.sum.ag.sample <- baypass.Mcali.sum.ag %>% 
    filter(SNP_id != unique(baypass.Mcali.sum.ag.out$chr)) %>% 
    filter(mean_M_P > baypass.g27343.mean.MP-0.15 & mean_M_P < baypass.g27343.mean.MP+0.15) |> 
    slice_sample(n = 1000)

# ================================================================================== #

# Correlation for sample SNPs

# Extract per site alpha for each of the 1,000 randomly sampled SNPs
Mcali.sample <- baypass.Mcali.sum %>% filter(SNP_id %in% baypass.Mcali.sum.ag.sample$SNP_id)

# Join with Mcali data
Mcali.sample <- left_join(Mcali.sample, Mcali, by="Site")
# Join with PC1 data
Mcali.sample <- left_join(Mcali.sample, PC1)

# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
baypass.cor.Mcali.sample.PC1.lat <- foreach(i=1:length(unique(Mcali.sample$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = Mcali.sample %>% filter(SNP_id == unique(Mcali.sample$SNP_id)[i])
    # Model
    mod.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + mean_integrated_thk, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(Mcali.sample$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "sample")
}
# Create column with absolute value of cor
baypass.cor.Mcali.sample.PC1.lat$abs.estimate <- abs(baypass.cor.Mcali.sample.PC1.lat$estimate)


# ================================================================================== #

# Graph correlation of outlier SNP compared to sample

# Set colors
cols <- c("#757474", "#29b3e6")

# Graph
pdf("output/figures/outlier_analyses/Correlation_abs_Mcali_baypass_PC1_lat_density_fill.pdf", width = 6, height = 6)
ggplot(baypass.cor.Mcali.sample.PC1.lat, aes(x = abs.estimate, fill = "group")) +
  geom_density(alpha = 0.7, lwd = 0.5) + 
  scale_fill_manual(values = cols[1]) + 
  geom_vline(xintercept = baypass.cor.Mcali.out.PC1.lat$abs.estimate, color = cols[2]) +
  xlim(0,1) + ylim(0,3) +
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()