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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)

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

# Load AFs for SNPs that beat POD in gene g27343 (produced in 05_graph_ph_outliers)
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
# ================================================================================== #

# Calculate AF for every SNP

# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# ================================================================================== #

# Extract and manipulate count and coverage

# Extract read count data for SNPs
ref_count <- pooldata@refallele.readcount
ref_count %>% as.data.frame -> count
names(count) = c(pooldata@poolnames)
count$SNP_id <- snp.info$SNP_id
count.melt <- reshape2::melt(count, id = "SNP_id", variable.name = "Site", value.name = "Count")

# Extract coverage data for SNPs
coverage <- pooldata@readcoverage
coverage %>% as.data.frame -> cov
names(cov) = c(pooldata@poolnames)
cov$SNP_id <- snp.info$SNP_id
cov.melt <- reshape2::melt(cov, id = "SNP_id", variable.name = "Site", value.name = "Cov")

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage
# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata@poolnames)
# Add SNP_id
afs$SNP_id <- snp.info$SNP_id
# Change format
afs.melt <- reshape2::melt(afs, id = "SNP_id", variable.name = "Site", value.name = "AF")

# Join
afs.all <- left_join(count.melt, cov.melt)
afs.all <- left_join(afs.all, afs.melt)
afs.all <- left_join(snp.info, afs.all)

# Save
save(afs.all, file = "data/processed/outlier_analyses/afs.all.RData")
load("data/processed/outlier_analyses/afs.all.RData")

# ================================================================================== #

# Summarize across sites to calculate global AFs
afs.all.ag <- afs.all %>% group_by(chr, pos, RefAllele, AltAllele, rs.id, SNP_id) %>%
    reframe(nSamps_poly=sum(AF>0 & AF<1, na.rm=T),
            nSamps_fixed=sum(AF==0 | AF==1, na.rm=T),
            global_af=sum(Count, na.rm=T)/sum(Cov, na.rm=T),
            mean_af=mean(AF, na.rm=T),
            poly_af=mean(AF[AF>0 & AF<1 & !is.na(AF)], na.rm=T))

# Save
save(afs.all.ag, file = "data/processed/outlier_analyses/afs.all.ag.RData")
load("data/processed/outlier_analyses/afs.all.ag.RData")

# ================================================================================== #

# Calculate correlation coef between AF and mean pH for 51 SNPs
#cor.g27343 <- foreach(i=1:length(unique(afs.ph.g27343.BF.POD$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = afs.ph.g27343.BF.POD %>% filter(SNP_id == unique(afs.ph.g27343.BF.POD$SNP_id)[i])
#    # Do correlation test
#    cor.i = cor.test(formula = ~ AF + ph_mean, data = tmp.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(afs.ph.g27343.BF.POD$SNP_id)[i],
#        statistic = cor.i[1],
#        p.value = cor.i[3],
#        estimate = cor.i[4],
#        group = "outliers")
#}
# Create column with absolute value of cor
#cor.g27343$abs.estimate <- abs(cor.g27343$estimate)


# Given that demography is important, account for that by using residuals from a regression with PC1 and AF
# Join with PC1 data
afs.ph.g27343.BF.POD <- left_join(afs.ph.g27343.BF.POD, PC1)
# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
tmp.cor.g27343.PC1.lat <- foreach(i=1:length(unique(afs.ph.g27343.BF.POD$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = afs.ph.g27343.BF.POD %>% filter(SNP_id == unique(afs.ph.g27343.BF.POD$SNP_id)[i])
    # Model
    mod.i = lm(AF~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + ph_mean, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(afs.ph.g27343.BF.POD$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "outliers")
}
# Create column with absolute value of cor
tmp.cor.g27343.PC1.lat$abs.estimate <- abs(tmp.cor.g27343.PC1.lat$estimate)

# ================================================================================== #

# Identify 1,000 SNPs with similar mean AF to those outlier SNPs in g27343

# Extract summary AF data for all outlier SNPs in g27343
afs.all.ag.g27343 <- afs.all.ag %>% filter(SNP_id %in% unique(afs.ph.g27343.BF.POD$SNP_id))
# Calc summary stats
g27343.mean.AF <- mean(afs.all.ag.g27343$mean_af)
g27343.sd.AF <- sd(afs.all.ag.g27343$mean_af)
g27343.range.AF <- range(afs.all.ag.g27343$mean_af)

# Extract a sample of 1000 SNPs that aren't in g27343 and have a mean AF within one sd of the mean AF for g27343
afs.all.ag.sample <- afs.all.ag %>% 
    filter(SNP_id != unique(afs.all.ag.g27343$chr)) %>% 
    filter(mean_af > g27343.mean.AF-g27343.sd.AF & mean_af < g27343.mean.AF+g27343.sd.AF) |> 
    slice_sample(n = 1000)

# ================================================================================== #

# Extract per site AF for each of the 1,000 randomly sampled SNPs
afs.sample <- afs.all %>% filter(SNP_id %in% afs.all.ag.sample$SNP_id)

# Join with pH data
afs.sample.ph <- left_join(afs.sample, ph, by="Site")

# Calculate correlation coef between AF and mean pH for sample of SNPs
#cor.sample <- foreach(i=1:length(unique(afs.sample.ph$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = afs.sample.ph %>% filter(SNP_id == unique(afs.sample.ph$SNP_id)[i])
#    # Do correlation test
#    cor.i = cor.test(formula = ~ AF + ph_mean, data = tmp.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(afs.sample.ph$SNP_id)[i],
#        statistic = cor.i[1],
#        p.value = cor.i[3],
#        estimate = cor.i[4],
#        group = "sample")
#}
# Create column with absolute value of cor
#cor.sample$abs.estimate <- abs(cor.sample$estimate)


# Given that demography is important, account for that by using residuals from a regression with PC1 and AF
# Join with PC1 data
afs.sample.ph <- left_join(afs.sample.ph, PC1)
# Calculate correlation coefficient between mean pH and the residuals from an AF model that accounted for PC1 (i.e., demography)
tmp.cor.sample.PC1.lat <- foreach(i=1:length(unique(afs.sample.ph$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = afs.sample.ph %>% filter(SNP_id == unique(afs.sample.ph$SNP_id)[i])
    # Model 
    mod.i = lm(AF~PC1+latitude, data = tmp.i)
    # Add to df
    tmp.i$residuals = mod.i$residuals
    # Do correlation test
    cor.i = cor.test(formula = ~ residuals + ph_mean, data = tmp.i)

    # Make table
    data.frame(
        SNP_id = unique(afs.sample.ph$SNP_id)[i],
        statistic = cor.i[1],
        p.value = cor.i[3],
        estimate = cor.i[4],
        group = "sample")
}
# Create column with absolute value of cor
tmp.cor.sample.PC1.lat$abs.estimate <- abs(tmp.cor.sample.PC1.lat$estimate)

# ================================================================================== #

# Join
#cor.join <- rbind(cor.g27343, cor.sample)

# Set colors
#cols <- c("orange", "#f53c3c")

# Graph Pearson's Correlation as density graph for outliers and sample
#pdf("output/figures/baypass/outliers/Correlation_ph_density.pdf", width = 12, height = 6)
#ggplot(cor.join, aes(x = estimate, colour = group)) +
#  geom_density(lwd = 1.2, linetype = 1) + 
#  scale_color_manual(values = cols) + 
#  labs(x = "Correlation", y = "Density") + 
#  theme_linedraw(base_size = 30) 
#dev.off()

#pdf("output/figures/baypass/outliers/Correlation_ph_density_fill.pdf", width = 12, height = 6)
#ggplot(cor.join, aes(x = estimate, fill = group)) +
#  geom_density(alpha = 0.7, color = NA) +
#  scale_fill_manual(values = cols) + 
#  labs(x = "Correlation", y = "Density") +
#  theme_linedraw(base_size = 30) 
#dev.off()

#pdf("output/figures/baypass/outliers/Correlation_abs_ph_density_fill.pdf", width = 12, height = 6)
#ggplot(cor.join, aes(x = abs.estimate, fill = group)) +
#  geom_density(alpha = 0.7, color = NA) +
#  scale_fill_manual(values = cols) + 
#  labs(x = "|Correlation|", y = "Density") +
#  theme_linedraw(base_size = 30) 
#dev.off()

# ================================================================================== #

# Join - PC1 incorporated
tmp.cor.PC1.join <- rbind(tmp.cor.g27343.PC1, tmp.cor.sample.PC1)

# Set colors
#cols <- c("#ff7b00", "#757474")
cols <- c("orange2", "#757474")

# Graph 
pdf("output/figures/baypass/outliers/Correlation_abs_ph_PC1_density_fill.pdf", width = 11, height = 4)
ggplot(tmp.cor.PC1.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,1) +
  scale_fill_manual(values = cols) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30)
dev.off()

# Join PC1 and lat
tmp.cor.PC1.lat.join <- rbind(tmp.cor.g27343.PC1.lat, tmp.cor.sample.PC1.lat)

# Graph
pdf("output/figures/baypass/outliers/Correlation_abs_ph_PC1_lat_density_fill.pdf", width = 11, height = 4)
ggplot(tmp.cor.PC1.lat.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,1) +
  scale_fill_manual(values = cols) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30)
dev.off()


# ================================================================================== #
# ================================================================================== #
# ================================================================================== #
# ================================================================================== #


# ================================================================================== #

# Extract data for just outlier SNPs on g27343

# NOTE baypass.ph.sum was created in src/05_outlier_analyses/06_corrected_xtx_ph.R
load("data/processed/baypass/abiotic/baypass.ph.sum.RData")

# Join outlier AF dataset with baypass corrected AF
ph.g27343 <- left_join(afs.ph.g27343.BF.POD, baypass.ph.sum, by = c("SNP_id", "Site"))
# Join with pH data
ph.g27343 <- left_join(ph.g27343, ph)

# Calculate correlation coef between AF and mean pH for 51 SNPs
#baypass.cor.g27343 <- foreach(i=1:length(unique(ph.g27343$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = ph.g27343 %>% filter(SNP_id == unique(ph.g27343$SNP_id)[i])
#    # Do correlation test
#    cor.i = cor.test(formula = ~ M_P.mean + ph_mean, data = tmp.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(ph.g27343$SNP_id)[i],
#        statistic = cor.i[1],
#        p.value = cor.i[3],
#        estimate = cor.i[4],
#        group = "outliers")
#}
# Create column with absolute value of cor
#baypass.cor.g27343$abs.estimate <- abs(baypass.cor.g27343$estimate)

#####

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

#####

# Extract per site AF for each of the 1,000 randomly sampled SNPs
ph.sample <- baypass.ph.sum %>% filter(SNP_id %in% baypass.ph.sum.ag.sample$SNP_id)

# Join with pH data
ph.sample <- left_join(ph.sample, ph, by="Site")

#####

# Calculate correlation coef between AF and mean pH for 51 SNPs
#baypass.cor.sample <- foreach(i=1:length(unique(ph.sample$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = ph.sample %>% filter(SNP_id == unique(ph.sample$SNP_id)[i])
#    # Do correlation test
#    cor.i = cor.test(formula = ~ M_P.mean + ph_mean, data = tmp.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(ph.sample$SNP_id)[i],
#        statistic = cor.i[1],
#        p.value = cor.i[3],
#        estimate = cor.i[4],
#        group = "sample")
#}
# Create column with absolute value of cor
#baypass.cor.sample$abs.estimate <- abs(baypass.cor.sample$estimate)

#####

# Join - PC1 incorporated
#baypass.cor.g27343.join <- rbind(baypass.cor.g27343, baypass.cor.sample)

# Set colors
cols <- c("#ff7b00", "#757474")

# Graph
#pdf("output/figures/baypass/outliers/Correlation_abs_ph_baypass_density_fill.pdf", width = 11, height = 4)
#ggplot(baypass.cor.g27343.join, aes(x = abs.estimate, fill = group)) +
#  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,1) +
#  scale_fill_manual(values = cols) + 
#  labs(x = "|Correlation|", y = "Density") +
#  theme_linedraw(base_size = 30)
#dev.off()


######

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

#####

# Join - PC1 and lat incorporated
baypass.cor.g27343.PC1.lat.join <- rbind(baypass.cor.g27343.PC1.lat, baypass.cor.sample.PC1.lat)

# Graph
pdf("output/figures/baypass/outliers/Correlation_abs_ph_baypass_PC1_lat_density_fill.pdf", width = 11, height = 4)
ggplot(baypass.cor.g27343.PC1.lat.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,0.5) +
  scale_fill_manual(values = cols) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30)
dev.off()

# Graph alt

cols_alt <- c(brewer.pal(7, "Oranges")[5], "#757474")
pdf("output/figures/baypass/outliers/Correlation_abs_ph_baypass_PC1_lat_density_fill_alt.pdf", width = 6, height = 6)
ggplot(baypass.cor.g27343.PC1.lat.join, aes(x = abs.estimate, fill = group)) +
  geom_density(alpha = 0.7, lwd = 0.5) + xlim(0,0.478) +
  scale_fill_manual(values = cols_alt) + 
  labs(x = "|Correlation|", y = "Density") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

