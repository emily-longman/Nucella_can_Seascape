# Redundancy analysis (RDA) as a genotype-environment association (GEA)
# Note: prior to running the R script, need to load R module on the VACC (module load R/4.4.1)
# Useful tutorial: https://popgen.nescent.org/2018-03-27_RDA_GEA.html

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
install.packages(c('poolfstat', 'data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'gameofthrones', 'vegan', 'psych'))
library(poolfstat)
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(vegan) # Used to run the RDA
library(psych)

# ================================================================================== #

# Load pooldata object
load("data/processed/pop_structure/pooldata.RData")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #
# ================================================================================== #

# Format read count data and subsample
# Note: RDA requires complete data frames (i.e., no missing data)

# Check structure
str(pooldata)

# Turn the ref allele read count data object into a dataframe and add SNP name (Chr and pos)
readcount <- as.data.frame(pooldata@refallele.readcount)
readcov <- as.data.frame(pooldata@readcoverage)
allele_freq <- as.data.frame(readcount/readcov)

# Extract SNP info
snp_name <- paste(pooldata@snp.info$Chromosome, pooldata@snp.info$Position)
rownames(allele_freq) <- snp_name

# Transpose dataframe
allele_freq_t <- t(allele_freq)

# Set rownames as sites
rownames(allele_freq_t) <- pooldata@poolnames

# Check to make sure no missing data
sum(is.na(allele_freq_t))

# Subsample SNP list for testing - 1,000,000 SNPs
allele_freq_t_sub <- allele_freq_t[,sample(1:ncol(allele_freq_t), 1000000)]

# ================================================================================== #

# Format bio-oracle data

# Check structure 
str(bio_oracle_sites_2010)

# Re-order bio-oracle data so in same order as read count data
sites <- rownames(allele_freq_t)
bio_oracle_sites_2010$location = factor(bio_oracle_sites_2010$location, levels = sites)
bio_oracle_sites_2010_ordered <- bio_oracle_sites_2010[order(bio_oracle_sites_2010$location), ]
# Make site names characters, not factors 
bio_oracle_sites_2010_ordered$location <- as.character(bio_oracle_sites_2010_ordered$location)

# Confirm that read count data and environmental data are in the same order
identical(rownames(allele_freq_t), bio_oracle_sites_2010_ordered[,13])

# ================================================================================== #

# Assess correlations among the Bio-oracle environmental variables
# Note: typically remove predictors with |r| > 0.7

# Bivariate scatter plots below the diagonal, histograms on the diagonal, and the Pearson correlation above the diagonal
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_correlations.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered[,4:12], scale=T)
dev.off()

# Many of the variables are correlated!!

# Do variable reduction. Remove one and then check, until remaining variables are not correlated

# Test removing variables
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_correlations_sub_test.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered[,c(7,8,9,10,11)], scale=T)
dev.off()

# One option is to keep: thetao_mean, chl_mean, O2_mean, and ph_min (that means cutting thetao_max, thetao_min, thetao_range, pH_mean, and so_mean)
bio_oracle_sites_2010_ordered_sub <- bio_oracle_sites_2010_ordered[,-c(1,2,3,4,5,6,11,13,14)]
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_correlations_sub.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered_sub[,1:5], scale=T)
dev.off()

# OR keep: thetao_mean, chl_mean, ph_min and ph_mean (that means cutting thetao_max, thetao_min, thetao_range, 02_mean, and so_mean; 
# Also ph_min and ph_mean are 0.71 cor
bio_oracle_sites_2010_ordered_sub2 <- bio_oracle_sites_2010_ordered[,-c(1,2,3,4,5,6,9,13,14)]
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_correlations_sub2.pdf", width = 10, height = 10)
pairs.panels(bio_oracle_sites_2010_ordered_sub2[,1:5], scale=T)
dev.off()

# ================================================================================== #
# ================================================================================== #

# Run the rda
N_can_rda <- rda(allele_freq_t ~ ., data = bio_oracle_sites_2010_ordered_sub, scale = T)
#N_can_rda_sub <- rda(readcount_t_sub ~ ., data = bio_oracle_sites_2010_ordered_sub, scale = T)

# ================================================================================== #

# Summary
N_can_rda
#Call: rda(formula = allele_freq_t ~ thetao_mean + chl_mean + o2_mean +
#ph_min + so_mean, data = bio_oracle_sites_2010_ordered_sub, scale = T)
#-- Model Summary --
#                Inertia Proportion Rank
#Total         8.277e+06  1.000e+00     
#Constrained   4.798e+06  5.796e-01    5
#Unconstrained 3.480e+06  4.204e-01   13
#Inertia is correlations
#-- Eigenvalues --
#Eigenvalues for constrained axes:
#   RDA1    RDA2    RDA3    RDA4    RDA5 
#3016365  676763  499643  364405  240386
#Eigenvalues for unconstrained axes:
#   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8    PC9   PC10   PC11 
#695221 483267 443629 337181 304638 236463 210065 188701 171075 158338 105783 
#  PC12   PC13 
# 79687  65594 

# Notes:
# Number of constrained RDA axes is the same as the number of predictors in the model.
# The proportion of the variance explained by the environmental predictors is given in the “Proportion” column of the model summary.

# Calculate R2 (must adjust based on the number of predictors)
RsquareAdj(N_can_rda)
# Our environmental predictors explain 41.79% of the variation! This is fairly high for the full SNP list since you would assume most are neutral.

# Calculate eigenvalues of constrained axes
summary(eigenvals(N_can_rda, model = "constrained"))

# Graph screeplot of canonical eigenvalues
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_scree.pdf", width = 10, height = 10)
screeplot(N_can_rda)
dev.off()
# The first constrained axes explains most of the variance by far!

# Test for significance of of SNP data and enviro predictors
sig_full <- anova.cca(N_can_rda, parallel=getOption("mc.cores")) # default is permutation=999
sig_full
# Takes a long time to finish

# Test for significance of each constrained axis (helps identify candidate loci)
signif.axis <- anova.cca(N_can_rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
# Cant get it to finish

# Check Variance Inflation Factors (VIF) - i.e., this tests for multicollinearity among predictors (you want them less than 10, ideally less than 5)
vif.cca(N_can_rda)
#thetao_mean    chl_mean     o2_mean      ph_min     so_mean 
#    2.090068    2.811177   11.244335    2.079426    8.231358 

# ================================================================================== #

# Graph the RDA
# Note: symmetrical scaling (i.e., scale the SNP and individual scores by the square root of the eigenvalues)

# SNPs are in red; the populations are the black circles; the blue vectors are the environmental predictors

# RDA1 vs RDA2
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_RDA1_RDA2.pdf", width = 10, height = 10)
plot(N_can_rda, scaling=3)
dev.off()

# RDA1 vs RDA3
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_RDA1_RDA3.pdf", width = 10, height = 10)
plot(N_can_rda, choices = c(1, 3), scaling=3)
dev.off()

# ================================================================================== #

# More detailed graphing

# Order levels/sites
levels(bio_oracle_sites_2010_ordered$location) <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
sites <- bio_oracle_sites_2010_ordered$location

# Set colors
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph RDA1 and RDA2
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_RDA1_RDA2_detailed.pdf", width = 10, height = 10)
plot(N_can_rda, type="n", scaling=3)
points(N_can_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
points(N_can_rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=mycolors)
text(N_can_rda, scaling=3, display="bp", col="#044ea7ea", cex=1)
legend("topleft", legend=levels(sites), bty="n", col="gray32", pch=21, cex=1, pt.bg=mycolors)
dev.off()

# Graph RDA1 and RDA3
pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_RDA1_RDA3_detailed.pdf", width = 10, height = 10)
plot(N_can_rda, type="n", scaling=3, choices=c(1,3))
points(N_can_rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
points(N_can_rda, display="sites", pch=21, cex=2, col="gray32", scaling=3, bg=mycolors)
text(N_can_rda, scaling=3, display="bp", col="#044ea7ea", cex=1)
legend("topleft", legend=levels(sites), bty="n", col="gray32", pch=21, cex=1, pt.bg=mycolors)
dev.off()

# ================================================================================== #

# Identify candidate SNPs

# SNP loadings ("species scores") for the first three constrained axes
load_rda <- scores(N_can_rda, choices=c(1:3), display="species")

pdf("output/figures/GEA/rda/Bio-oracle/Bio-oracle_RDA_loadings.pdf", width = 12, height = 6)
par(mfrow=c(1,3))
hist(load_rda[,1], main="Loadings on RDA1")
hist(load_rda[,2], main="Loadings on RDA2")
hist(load_rda[,3], main="Loadings on RDA3")
dev.off()

# SNPs loadings for RDA 1 are showing a relationship with the environmental predictors! Thus, it is not just a few outlier loci that are under selection
# This is true for both subsets of enviro predictors

# ================================================================================== #
# ================================================================================== #

# Likely shouldn't proceed, but just curious to see

# Outlier function (find loadinds that are +/- z sd from the mean loading)
outliers <- function(x,z){
lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/- z sd from mean loading     
x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

# Identify outliers for each rda axis
cand1 <- outliers(load_rda[,1],3.5) #0
cand2 <- outliers(load_rda[,2],3.5) #494
cand3 <- outliers(load_rda[,3],3.5) #481

# Total number of outliers
ncand <- length(cand1) + length(cand2) + length(cand3)
ncand #975

# Format data (dataframe with axis, SNP name, loadinga nd correlation with predictor)
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)

# Add correlations of each candidate SNP with the enviro predictors
df <- matrix(nrow=(ncand), ncol=5)  # 5 columns for 5 predictors
colnames(df) <- c("thetao_mean", "chl_mean", "o2_mean", "ph_min", "so_mean")

# Problems with this bc snps names don't match - for cand$snp they are based on the column number e.g., "col13358"
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- readcount_t[,nam]
  df[i,] <- apply(bio_oracle_sites_2010_ordered_sub, 2 ,function(x) cor(x, snp.gen))
}

# Combine data frame
cand <- cbind.data.frame(cand, df)
head(cand)