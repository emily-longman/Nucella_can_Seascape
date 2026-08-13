# Graph Mcali outliers

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
bf.McaliIntThk.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", header=T)

# Extract those with BF>15
bf.Mcali.sub <- bf.McaliIntThk.mean.sum.outliers.annotated %>% filter(bf_db.mean >= 15)

# Load Mcali data
Mcalifornianus_data_clean <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean.csv", header=T)

# Extract just mean pH and rename location to Site
Mcalifornianus_data_clean <- Mcalifornianus_data_clean[,c(1,3,4,7)]
Mcalifornianus_data_clean <- Mcalifornianus_data_clean %>% rename(Site = Site.Code)

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
selected_SNPs_Mcali <- snp.info %>% filter(SNP_id %in% bf.Mcali.sub$SNP_id)
# Get index of SNPs
selected_SNPs_Mcali <- as.integer(selected_SNPs_Mcali$rs.id)
# Subset the pooldata object using the selected SNP indices
pooldata_Mcali <- pooldata.subset(pooldata, snp.index = selected_SNPs_Mcali)

# Extract and manipulate snp info for significant SNPs

# Extract SNP info for significant SNPs
pooldata_Mcali@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info.Mcali
# Rename columns
names(snp.info.Mcali)[1:2] = c("chr","pos")
# Make snp_id column
snp.info.Mcali %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info.Mcali

# ================================================================================== #

# Extract and manipulate coverage for significant SNPs

# Extract read count and coverage data for SNPs
ref_count <- pooldata_Mcali@refallele.readcount
coverage <- pooldata_Mcali@readcoverage

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage

# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata_Mcali@poolnames)

# Add SNP_id
afs$SNP_id <- snp.info.Mcali$SNP_id

# Change format
afs.melt <- reshape2::melt(afs, id = "SNP_id", variable.name = "Site", value.name = "AF")

# Join with pH data
afs.Mcali <- left_join(afs.melt, Mcalifornianus_data_clean, by="Site")

# Remove rows w/ NAs - since not all sites have drilling data
afs.Mcali <- na.omit(afs.Mcali)

# Make Site an ordered factor
lat.order <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
afs.Mcali <- afs.Mcali %>% mutate(Site = factor(Site, levels = lat.order))

# Color palette
viridiscolors <- viridis(n=19)

# Subset so only SNPs with BF > 20 or >18
afs.Mcali.BF <- afs.Mcali %>% filter(SNP_id %in% bf.Mcali.sub[which(bf.Mcali.sub$bf_db.mean > 23),]$SNP_id)

# Graph and color by Site
pdf("output/figures/baypass/outliers/Mcali_BF23.pdf", width = 18, height = 12)
ggplot(afs.Mcali.BF, aes(x=AF, y=mean_integrated_thk, fill=Site)) +
  geom_point(alpha=0.6, size = 4, shape = 21) + labs(x="Allele Frequency", y="Mcali Thickness") +
  facet_wrap(~SNP_id) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 26) + theme(legend.position="none")
dev.off()
