# Compare effect sizes (Cohen's F2) for outlier SNPs associated with mean pH and those associated with M. californianus cross sectional thickness with a random sample of 1,000 SNPs

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

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/outliers")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load AF per SNP
load("data/processed/outlier_analyses/afs.all.RData")

# pH Data
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
# M. cali thk
load("data/processed/baypass/biotic/bf.Mcali.IntThk.sum.Rdata")

# Load thresholds
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")
bf.POD.thr.ph <- bf.POD.thr
load("data/processed/baypass/biotic/Mcali_IntegratedThk_POD_thr.Rdata")
bf.POD.thr.McaliThk <- bf.POD.thr

# Load Baypass alpha(ij)
load("data/processed/baypass/abiotic/baypass.ph.sum.RData")
load("data/processed/baypass/biotic/baypass.Mcali.sum.RData")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)
# Extract just mean pH and rename location to Site
ph <- bio_oracle_sites_2010[,c(1,2,3,11,13)]
ph <- ph %>% rename(Site = location)

# Load M californianus shell thickness data
Mcalifornianus_data <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean_18pop.csv", header=T)
Mcali <- Mcalifornianus_data[,c(1,3,4,7)]
Mcali <- Mcali %>% rename(Site = Site.Code, latitude = Latitude, longitude = Longitude)

# Load PCA data for demography
pca.df <- read.csv("data/processed/outlier_analyses/pca.csv")
colnames(pca.df)[1] <- "Site"
PC1 <- pca.df[,c(1,2)]

# ================================================================================== #

# For both abiotic and biotic extract SNPs that beat threshold

# Mean pH
bf.ph.mean.sum.outliers <- bf.ph.mean.sum %>% 
    filter(bf_db.mean > bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)]) %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# M cali thk
bf.McaliIntThk.mean.sum.outliers <- bf.McaliIntThk.mean.sum %>% 
    filter(bf_db.mean > bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)]) %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Extract AF data for each SNP subset

# mean pH
ph.outliers.af <- afs.all %>% filter(SNP_id %in% bf.ph.mean.sum.outliers$SNP_id)
ph.outliers.af <- left_join(ph.outliers.af, ph)

# M cali thickness (also drop ARA site)
Mcali.outliers <- afs.all %>% filter(SNP_id %in% bf.McaliIntThk.mean.sum.outliers$SNP_id)
Mcali.outliers <- left_join(Mcali.outliers, Mcali) %>% drop_na(mean_integrated_thk)

# ================================================================================== #

# Calculate effect size (Cohen's F squared)

# Abiotic

# Join outlier SNPs with baypass corrected AFs data
ph.outliers.baypass <- left_join(ph.outliers.af, baypass.ph.sum, by = join_by(chr, pos, SNP_id, Site)) 
ph.outliers.baypass <- left_join(ph.outliers.baypass, PC1)
# Calculate local effect size (Cohen's f2)
ph.cohens.fsq.PC1.lat <- foreach(i=1:length(unique(ph.outliers.baypass$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = ph.outliers.baypass %>% filter(SNP_id == unique(ph.outliers.baypass$SNP_id)[i])
    # Model - full
    mod.full.i = lm(M_P.mean~ph_mean+PC1+latitude, data = tmp.i)
    # Extract r2 - full mod
    r2.full.i = summary(mod.full.i)$r.squared
    # Model - reduced
    mod.reduced.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Extract r2 - full mod
    r2.reduced.i = summary(mod.reduced.i)$r.squared
    # Calculate effect size (Cohen's F sq)
    fsq.i = (r2.full.i-r2.reduced.i)/(1-r2.full.i)

    # Make table
    data.frame(
        SNP_id = unique(ph.outliers.baypass$SNP_id)[i],
        cohens.fsq = fsq.i,
        group = "ph")
}

#######

# Biotic

# Join outlier SNPs with baypass corrected AFs data
Mcali.outliers.baypass <- left_join(Mcali.outliers, baypass.Mcali.sum, by = join_by(chr, pos, SNP_id, Site)) 
Mcali.outliers.baypass <- left_join(Mcali.outliers.baypass, PC1)
# Calculate local effect size (Cohen's f2)
Mcali.cohens.fsq.PC1.lat <- foreach(i=1:length(unique(Mcali.outliers.baypass$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = Mcali.outliers.baypass %>% filter(SNP_id == unique(Mcali.outliers.baypass$SNP_id)[i])
    # Model - full
    mod.full.i = lm(M_P.mean~mean_integrated_thk+PC1+latitude, data = tmp.i)
    # Extract r2 - full mod
    r2.full.i = summary(mod.full.i)$r.squared
    # Model - reduced
    mod.reduced.i = lm(M_P.mean~PC1+latitude, data = tmp.i)
    # Extract r2 - full mod
    r2.reduced.i = summary(mod.reduced.i)$r.squared
    # Calculate effect size (Cohen's F sq)
    fsq.i = (r2.full.i-r2.reduced.i)/(1-r2.full.i)

    # Make table
    data.frame(
        SNP_id = unique(Mcali.outliers$SNP_id)[i],
        cohens.fsq = fsq.i,
        group = "Mcali")
}

# ================================================================================== #

# Join with BF

# Join with BF - alpha
ph.cohens.fsq.PC1.lat.bf <- left_join(bf.ph.mean.sum.outliers, ph.cohens.fsq.PC1.lat)
Mcali.cohens.fsq.PC1.lat.bf <- left_join(bf.McaliIntThk.mean.sum.outliers, Mcali.cohens.fsq.PC1.lat)

# Save M cali data for following script
save(Mcali.cohens.fsq.PC1.lat.bf, file = "data/processed/baypass/biotic/Mcali.cohens.fsq.PC1.lat.bf.RData")

# ================================================================================== #

# Graph - BF with cohen's F sq

# Abiotic

# Graph BF vs Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_BF_ph_baypass_PC1_lat.pdf", width = 8, height = 6)
ggplot(ph.cohens.fsq.PC1.lat.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.6) + 
  xlim(bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], 30.5) +
  labs(x = "BF mean pH", y = "Cohen's f2") + 
  theme_linedraw(base_size = 30)
dev.off()
# Graph and order top SNPs by Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_ph_order_baypass_PC1_lat.pdf", width = 8, height = 6)
ggplot(ph.cohens.fsq.PC1.lat.bf, aes(x = reorder(SNP_id, cohens.fsq), y = cohens.fsq)) +
  geom_point() + ylim(0,17.5) +
  labs(y = "Cohen's f2", x = "Top pH SNPs") + 
  theme_linedraw(base_size = 30) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()

# Biotic

# Set colors
#cols <- c("#757474", "#29b3e6")
cols <- c("#757474", "#9e03a9")

# Graphs with bypass alpha
# Graph BF vs Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_BF_Mcali_baypass_PC1_lat.pdf", width = 6.3, height = 6)
ggplot(Mcali.cohens.fsq.PC1.lat.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.7, size = 6, aes(color = cohens.fsq >= 8 & bf_db.mean >= 25)) + 
  scale_color_manual(values = cols) +
  ylim(0,16.5) + xlim(bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], 30.5) +
  labs(x = "BF", y = expression(paste("Cohen's ", italic("f"^2)))) + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/outlier_analyses/Cohensf2_log2_BF_Mcali_baypass_PC1_lat.pdf", width = 8, height = 6)
ggplot(Mcali.cohens.fsq.PC1.lat.bf, aes(x = bf_db.mean, y = log2(cohens.fsq))) +
  geom_point(alpha = 0.6, size = 6, aes(color = cohens.fsq >= 8 & bf_db.mean >= 25)) + 
  scale_color_manual(values = cols) +
  geom_hline(yintercept = log2(0.02), col="red") +
  geom_hline(yintercept = log2(0.2), col="red") +
  geom_hline(yintercept = log2(0.4), col="red") +
  #ylim(0,16.5) + xlim(bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], 30.5) +
  labs(x = "BF", y = "Log2(Cohen's f2)") + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph and order top SNPs by Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_Mcali_order_baypass_PC1_lat.pdf", width = 8, height = 6)
ggplot(Mcali.cohens.fsq.PC1.lat.bf, aes(x = reorder(SNP_id, cohens.fsq), y = cohens.fsq)) +
  geom_point() + ylim(0,17.5) +
  labs(y = "Cohen's f2", x = "Top McaliThk SNPs") + 
  theme_linedraw(base_size = 30) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()


# ================================================================================== #

# Extract outlier (ntLink_3488_113935)
#Mcali.top.outlier <- Mcali.cohens.fsq.PC1.bf[which(Mcali.cohens.fsq.PC1.bf$cohens.fsq > 8),]
Mcali.top.outlier <- Mcali.cohens.fsq.PC1.lat.bf[which(Mcali.cohens.fsq.PC1.lat.bf$cohens.fsq > 8 & Mcali.cohens.fsq.PC1.lat.bf$bf_db.mean >= 25),]
#Mcali.cohensf2.outliers <- Mcali.cohens.fsq.PC1.lat.bf[which(Mcali.cohens.fsq.PC1.lat.bf$cohens.fsq > 11),]

# Get AF for top outlier SNP
Mcali.outliers.outlier <- Mcali.outliers[which(Mcali.outliers$SNP_id == Mcali.top.outlier$SNP_id),]
write.csv(Mcali.outliers.outlier, file = "data/processed/baypass/afs.McaliThk.outlier.csv", row.names=F)

# Make Site an ordered factor
lat.order <- c("FC", "SLR", "SH", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
Mcali.outliers.outlier <- Mcali.outliers.outlier %>% mutate(Site = factor(Site, levels = lat.order))
# Add column that highlights N, and S
Mcali.outliers.outlier <- Mcali.outliers.outlier %>% mutate(shape = case_when(Site %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   Site %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "SH", "SLR", "FC") ~ "N"))

# Color palette
viridiscolors <- viridis(n=19)
viridiscolors <- viridiscolors[-4]

# Graph AF vs Mcali thickness
pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk.pdf", width = 6.28, height = 6)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = Site)) +
  geom_point(alpha=0.8, size = 9) + scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  labs(x="Allele Frequency", y="Thickness") + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

Mcali.outliers.outlier <- left_join(Mcali.outliers.outlier, PC1)
pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk_PCA.pdf", width = 6.28, height = 6)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = PC1)) +
  geom_point(alpha=0.8, size = 9) + scale_shape_manual(values = c(21, 23)) + #scale_fill_manual(values = viridiscolors) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  scale_fill_gradientn(colours=brewer.pal(9, "RdGy")) +
  labs(x="Allele Frequency", y="Thickness") + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk_legend.pdf", width = 9, height = 9)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = Site)) +
  geom_point(alpha=0.8, size = 6) + scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) + 
  labs(x="Allele Frequency", y=expression(paste(italic("M. californianus"), " thickness"))) + 
  theme_linedraw(base_size = 30) + 
  guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none")
dev.off()


pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk_legend.pdf", width = 9, height = 9)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = Site)) +
  geom_point(alpha=0.8, size = 6) + scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) + 
  labs(x="Allele Frequency", y=expression(paste(italic("M. californianus"), " thickness"))) + 
  theme_linedraw(base_size = 30) + 
  guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none")
dev.off()

