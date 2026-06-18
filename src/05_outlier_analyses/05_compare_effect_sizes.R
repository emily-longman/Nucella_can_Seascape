# Compare effect sizes for outlier SNPs associated with mean pH and those associated with M. californianus cross sectional thickness

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis', 'stats', 'fastglm'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)
library(stats)
library(fastglm)

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
#ph.cohens.fsq <- foreach(i=1:length(unique(ph.outliers.af$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = ph.outliers.af %>% filter(SNP_id == unique(ph.outliers.af$SNP_id)[i])
#    # Model
#    mod.i = lm(AF~ph_mean, data = tmp.i)
#    # Extract r2
#    r2.i = summary(mod.i)$r.squared
#    # Calculate effect size (Cohen's F sq)
#    fsq.i = r2.i/(1-r2.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(ph.outliers.af$SNP_id)[i],
#        cohens.fsq = fsq.i,
#        group = "ph")
#}


ph.outliers.af <- left_join(ph.outliers.af, PC1)
ph.cohens.fsq.PC1 <- foreach(i=1:length(unique(ph.outliers.af$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = ph.outliers.af %>% filter(SNP_id == unique(ph.outliers.af$SNP_id)[i])
    # Model - full
    mod.full.i = lm(AF~ph_mean+PC1, data = tmp.i)
    # Extract r2 - full mod
    r2.full.i = summary(mod.full.i)$r.squared
    # Model - reduced
    mod.reduced.i = lm(AF~PC1, data = tmp.i)
    # Extract r2 - full mod
    r2.reduced.i = summary(mod.reduced.i)$r.squared
    # Calculate effect size (Cohen's F sq)
    fsq.i = (r2.full.i-r2.reduced.i)/(1-r2.full.i)

    # Make table
    data.frame(
        SNP_id = unique(ph.outliers.af$SNP_id)[i],
        cohens.fsq = fsq.i,
        group = "ph")
}


# Use baypass corrected AF
# Extract outliers and join with baypass data
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

# Biotic
#Mcali.cohens.fsq <- foreach(i=1:length(unique(Mcali.outliers$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
#    
#    # Extract rows associated with SNP
#    tmp.i = Mcali.outliers %>% filter(SNP_id == unique(Mcali.outliers$SNP_id)[i])
#    # Model
#    mod.i = lm(AF~mean_integrated_thk, data = tmp.i)
#    # Extract r2
#    r2.i = summary(mod.i)$r.squared
#    # Calculate effect size (Cohen's F sq)
#    fsq.i = r2.i/(1-r2.i)
#
#    # Make table
#    data.frame(
#        SNP_id = unique(Mcali.outliers$SNP_id)[i],
#        cohens.fsq = fsq.i,
#        group = "Mcali")
#}


Mcali.outliers <- left_join(Mcali.outliers, PC1)
Mcali.cohens.fsq.PC1 <- foreach(i=1:length(unique(Mcali.outliers$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = Mcali.outliers %>% filter(SNP_id == unique(Mcali.outliers$SNP_id)[i])
    # Model - full
    mod.full.i = lm(AF~mean_integrated_thk+PC1, data = tmp.i)
    # Extract r2 - full mod
    r2.full.i = summary(mod.full.i)$r.squared
    # Model - reduced
    mod.reduced.i = lm(AF~PC1, data = tmp.i)
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


# Use baypass corrected AF
# Extract outliers and join with baypass data
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
#ph.cohens.fsq.bf <- left_join(bf.ph.mean.sum.outliers, ph.cohens.fsq)
#Mcali.cohens.fsq.bf <- left_join(bf.McaliIntThk.mean.sum.outliers, Mcali.cohens.fsq)

ph.cohens.fsq.PC1.bf <- left_join(bf.ph.mean.sum.outliers, ph.cohens.fsq.PC1)
Mcali.cohens.fsq.PC1.bf <- left_join(bf.McaliIntThk.mean.sum.outliers, Mcali.cohens.fsq.PC1)

ph.cohens.fsq.PC1.lat.bf <- left_join(bf.ph.mean.sum.outliers, ph.cohens.fsq.PC1.lat)
Mcali.cohens.fsq.PC1.lat.bf <- left_join(bf.McaliIntThk.mean.sum.outliers, Mcali.cohens.fsq.PC1.lat)

# ================================================================================== #

# Graph - BF with cohen's F sq

# Abiotic
#pdf("output/figures/outlier_analyses/CohensF_BF_ph.pdf", width = 8, height = 6)
#ggplot(ph.cohens.fsq.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
#  geom_point(alpha = 0.6) + 
#  ylim(0,10) + xlim(bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], 30.5) +
#  labs(x = "pH BF", y = "Cohen's f2") + 
#  theme_linedraw(base_size = 30)
#dev.off()

# Graph BF vs Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_BF_ph_PC1.pdf", width = 8, height = 6)
ggplot(ph.cohens.fsq.PC1.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.6) + 
  xlim(bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], 30.5) +
  labs(x = "BF mean pH", y = "Cohen's f2") + 
  theme_linedraw(base_size = 30)
dev.off()

# Graph and order top SNPs by Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_ph_order_PC1.pdf", width = 8, height = 6)
ggplot(ph.cohens.fsq.PC1, aes(x = reorder(SNP_id, cohens.fsq), y = cohens.fsq)) +
  geom_point() + ylim(0,24) +
  labs(y = "Cohen's f2", x = "Top pH SNPs") + 
  theme_linedraw(base_size = 30) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()


# Graphs with bypass alpha
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
#pdf("output/figures/outlier_analyses/CohensF_BF_Mcali.pdf", width = 8, height = 6)
#ggplot(Mcali.cohens.fsq.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
#  geom_point(alpha = 0.6) + 
#  ylim(0,10) + xlim(bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], 30.5) +
#  labs(x = "Mcali BF", y = "Cohen's F2") + 
#  theme_linedraw(base_size = 30)
#dev.off()

# Set colors
cols <- c("#757474", "#29b3e6")

# Graph BF vs Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_BF_Mcali_PC1.pdf", width = 8, height = 6)
ggplot(Mcali.cohens.fsq.PC1.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.6, size = 6, aes(color = cohens.fsq >= 8)) + 
  scale_color_manual(values = cols) +
  ylim(0,10) + xlim(bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], 30.5) +
  labs(x = "BF", y = "Cohen's f2") + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph and order top SNPs by Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_Mcali_order_PC1.pdf", width = 8, height = 6)
ggplot(Mcali.cohens.fsq.PC1, aes(x = reorder(SNP_id, cohens.fsq), y = cohens.fsq)) +
  geom_point() + ylim(0,24) +
  labs(y = "Cohen's f2", x = "Top McaliThk SNPs") + 
  theme_linedraw(base_size = 30) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()

# Graphs with bypass alpha
# Graph BF vs Cohen's f2
pdf("output/figures/outlier_analyses/Cohensf2_BF_Mcali_baypass_PC1_lat.pdf", width = 8, height = 6)
ggplot(Mcali.cohens.fsq.PC1.lat.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.6, size = 6, aes(color = cohens.fsq >= 8 & bf_db.mean >= 25)) + 
  scale_color_manual(values = cols) +
  ylim(0,16.5) + xlim(bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], 30.5) +
  labs(x = "BF", y = "Cohen's f2") + 
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
Mcali.top.outlier <- Mcali.cohens.fsq.PC1.bf[which(Mcali.cohens.fsq.PC1.bf$cohens.fsq > 8),]
Mcali.top.outlier <- Mcali.cohens.fsq.PC1.lat.bf[which(Mcali.cohens.fsq.PC1.lat.bf$cohens.fsq > 8 & Mcali.cohens.fsq.PC1.lat.bf$bf_db.mean >= 25),]
Mcali.cohensf2.outliers <- Mcali.cohens.fsq.PC1.lat.bf[which(Mcali.cohens.fsq.PC1.lat.bf$cohens.fsq > 11),]

# Get AF for top outlier SNP
Mcali.outliers.outlier <- Mcali.outliers[which(Mcali.outliers$SNP_id == Mcali.top.outlier$SNP_id),]


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
pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk.pdf", width = 8, height = 6)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = Site)) +
  geom_point(alpha=0.8, size = 9) + scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) + 
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  labs(x="Allele Frequency", y="Thickness") + 
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/outlier_analyses/Mcali_AF_McaliThk_legend.pdf", width = 9, height = 7.5)
ggplot(Mcali.outliers.outlier, aes(x = AF, y = mean_integrated_thk, shape = shape, fill = Site)) +
  geom_point(alpha=0.8, size = 6) + scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) + 
  labs(x="Allele Frequency", y=expression(paste(italic("M. californianus"), " thickness"))) + 
  theme_linedraw(base_size = 30) + 
  guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none")
dev.off()



# ================================================================================== #
# ================================================================================== #
# ================================================================================== #
# ================================================================================== #






# ================================================================================== #

######## TEST GLM
# Doesn't make the most sense bc doesn't have a clear r2 so can't calculate an effect size easily
# Define anova function for comparing glms
anovaFun <- function(m1, m2) {
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  chisq <- -2*(ll1-ll2)
  1-pchisq(chisq, parameter)
    }

# GlM approach to include nEff in model
ph.glm <- foreach(i=1:length(unique(ph.outliers.af$SNP_id)), .combine = "rbind", .errorhandling = "remove")%do%{
    
    # Extract rows associated with SNP
    tmp.i = ph.outliers.af %>% filter(SNP_id == unique(ph.outliers.af$SNP_id)[i])

    # Sample size of each pool
    nSnail=20
    # Calculate the mean effective coverage ('nEff') (note: each pool consists of 20 dogwhelks)
    tmp.i <- tmp.i %>% mutate(nEff = round((Cov*2*nSnail)/(2*nSnail+Cov-1)))
    # Calculate the effective allele freq
    tmp.i <- tmp.i %>% mutate(AF_nEff = round(AF*nEff)/nEff)
    # Model allele freq
    # GLMs - a reduced model with just demography (PC1) and a model with demography and mean pH
    y <- tmp.i$AF_nEff
    X.reduced <- model.matrix(~PC1, tmp.i)
    X.full <- model.matrix(~PC1+ph_mean, tmp.i)
    t1.reduced <- fastglm(x=X.reduced, y=y, family=binomial(), weights=tmp.i$nEff, method=0)
    t1.full <- fastglm(x=X.full, y=y, family=binomial(), weights=tmp.i$nEff, method=0)

    # Make table
    data.frame(
        SNP_id = unique(ph.outliers.af$SNP_id)[i],
        AIC_reduced = c(AIC(t1.reduced)),
        AIC_full = c(AIC(t1.full)),
        p_lrt = anovaFun(t1.reduced, t1.full),
        group = "ph")
}
ph.glm <- ph.glm %>% mutate(AIC_delta = AIC_reduced-AIC_full)

ph.glm.bf <- left_join(bf.ph.mean.sum.outliers, ph.glm)

# Graph BF vs delta AIC
pdf("output/figures/outlier_analyses/Delta_AIC_BF_ph.pdf", width = 8, height = 6)
ggplot(ph.glm.bf, aes(x = bf_db.mean, y = cohens.fsq)) +
  geom_point(alpha = 0.6) + 
  xlim(bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], 30.5) +
  labs(x = "BF mean pH", y = "Cohen's f2") + 
  theme_linedraw(base_size = 30)
dev.off()

# Graph and order top SNPs by delta AIC
pdf("output/figures/outlier_analyses/Delta_ph_order.pdf", width = 8, height = 6)
ggplot(ph.cohens.fsq.PC1, aes(x = reorder(SNP_id, cohens.fsq), y = cohens.fsq)) +
  geom_point() + ylim(0,24) +
  labs(y = "Cohen's f2", x = "Top pH SNPs") + 
  theme_linedraw(base_size = 30) + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank())
dev.off()

######## TEST GLM