# Graph pH outliers

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
out_fig_dir <- paste("output/figures/baypass/outliers")
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

# Remove rows w/ NAs - since not all sites have drilling data
afs.ph <- na.omit(afs.ph)

# Make Site an ordered factor
lat.order <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
afs.ph <- afs.ph %>% mutate(Site = factor(Site, levels = lat.order))

# Color palette
viridiscolors <- viridis(n=19)

# Subset so only SNPs with BF > 20 or >18
afs.ph.BF <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 20),]$SNP_id)
afs.ph.BF19 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 19),]$SNP_id)
afs.ph.BF18 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.g27343[which(bf.ph.mean.g27343$bf_db.mean > 18),]$SNP_id)

# Graph and color by Site
pdf("output/figures/baypass/outliers/pH_g27343_BF20.pdf", width = 18, height = 4)
ggplot(afs.ph.BF, aes(x=AF, y=ph_mean, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency", y="pH") +
  facet_wrap(~SNP_id, ncol = 4) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 26) + theme(legend.position="none")
dev.off()

# Graph and color by Site
pdf("output/figures/baypass/outliers/pH_g27343_BF19.pdf", width = 14, height = 9)
ggplot(afs.ph.BF19, aes(x=AF, y=ph_mean, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Freq", y="pH") +
  facet_wrap(~SNP_id) + scale_fill_manual(values = viridiscolors) +
  theme_bw(base_size = 26)
dev.off()

# Graph and color by Site
pdf("output/figures/baypass/outliers/pH_g27343_BF18.pdf", width = 16, height = 12)
ggplot(afs.ph.BF18, aes(x=AF, y=ph_mean, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Freq", y="pH") +
  facet_wrap(~SNP_id) + scale_fill_manual(values = viridiscolors) +
  theme_bw(base_size = 26)
dev.off()
