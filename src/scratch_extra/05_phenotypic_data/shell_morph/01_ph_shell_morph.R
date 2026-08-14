# Analyse top ph SNPs with morphology

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/phenotypic_data")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load SNPs of interest
bf.ph.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header=T)
# Extract just those on g27343
bf.ph.mean.g27343 <- bf.ph.mean.sum.outliers.annotated %>% filter(Gene_Name == "g27343")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Extract just mean pH and rename location to Site
ph <- bio_oracle_sites_2010[,c(1,2,3,11,13)]
ph <- ph %>% rename(Site = location)

# Morphology data
morph <- fread("data/raw/phenotypic_data/pc.morphology.csv")
morph <- morph %>% rename(Site = Site.Code)

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Subset pooldata for Baypass outlier SNPs
# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Filter pooldata for Baypass SNPs of interest
selected_SNPs_ph <- snp.info %>% filter(SNP_id %in% bf.ph.mean.g27343$SNP_id)
# Get index of SNPs
selected_SNPs_ph_index <- as.integer(selected_SNPs_ph$rs.id)
# Subset the pooldata object using the selected SNP indices
pooldata_ph <- pooldata.subset(pooldata, snp.index = selected_SNPs_ph_index)

# Extract and manipulate snp info for significant SNPs

# Extract SNP info for significant SNPs
pooldata_ph@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info.ph
# Rename columns
names(snp.info.ph)[1:2] = c("chr","pos")
# Make snp_id column
snp.info.ph %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info.ph

# ================================================================================== #

# Extract and manipulate coverage for significant SNPs

# Extract read count and coverage data for SNPs
ref_count <- pooldata_ph@refallele.readcount
coverage <- pooldata_ph@readcoverage

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage

# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata_ph@poolnames)

# Add SNP_id
afs$SNP_id <- snp.info.ph$SNP_id

# Change format
afs.melt <- reshape2::melt(afs, id = "SNP_id", variable.name = "Site", value.name = "AF")

# Join with pH data
afs.ph <- left_join(afs.melt, ph, by="Site")
# Join with morph data
afs.ph.morph <- left_join(afs.ph, morph, by="Site")

# Remove rows w/ NAs - since not all sites
afs.ph.morph <- na.omit(afs.ph.morph)

# Make Site an ordered factor
lat.order <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
afs.ph.morph <- afs.ph.morph %>% mutate(Site = factor(Site, levels = lat.order))

# Color palette
viridiscolors <- viridis(n=19)

# Subset so only SNPs with BF > 20 or >18
afs.ph.morph.BF.20 <- afs.ph.morph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 20),]$SNP_id)
afs.ph.morph.BF19 <- afs.ph.morph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 19),]$SNP_id)
afs.ph.morph.BF18 <- afs.ph.morph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 18),]$SNP_id)


# Graph AF vs morph - mean PC 1
pdf("output/figures/phenotypic_data/pH_g27343_BF19_shell_morph_mean_pc1.pdf", width = 12, height = 10)
ggplot(afs.ph.morph.BF19, aes(x=AF, y=mean.pc1, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency") +
  facet_wrap(~SNP_id, ncol = 3) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 30) + theme(legend.position="none")
dev.off()

# Graph AF vs morph - mean PC 2
pdf("output/figures/phenotypic_data/pH_g27343_BF19_shell_morph_mean_pc2.pdf", width = 12, height = 10)
ggplot(afs.ph.morph.BF19, aes(x=AF, y=mean.pc2, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency") +
  facet_wrap(~SNP_id, ncol = 3) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 30) + theme(legend.position="none")
dev.off()

# Graph AF vs morph - CV PC 1 
pdf("output/figures/phenotypic_data/pH_g27343_BF19_shell_morph_CV_pc1.pdf", width = 12, height = 10)
ggplot(afs.ph.morph.BF19, aes(x=AF, y=CV.pc1, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency") +
  facet_wrap(~SNP_id, ncol = 3) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 30) + theme(legend.position="none")
dev.off()

# Graph AF vs morph - CV PC 2
pdf("output/figures/phenotypic_data/pH_g27343_BF19_shell_morph_CV_pc2.pdf", width = 12, height = 10)
ggplot(afs.ph.morph.BF19, aes(x=AF, y=CV.pc2, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency") +
  facet_wrap(~SNP_id, ncol = 3) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 30) + theme(legend.position="none")
dev.off()


####

# Graph mean PC1 and mean PC2 with AF for top SNPs
pdf("output/figures/phenotypic_data/pH_g27343_BF19_shell_morph_AF.pdf", width = 14, height = 10)
ggplot(afs.ph.morph.BF19, aes(x=mean.pc1, y=mean.pc2, fill=AF)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + 
  facet_wrap(~SNP_id, ncol = 3) + 
  geom_hline(yintercept = 0) + 
  geom_vline(xintercept = 0) + 
  scale_fill_gradient2(low = "steelblue", ,high = "firebrick", midpoint = 0.5) +
  theme_bw(base_size = 30)
dev.off()