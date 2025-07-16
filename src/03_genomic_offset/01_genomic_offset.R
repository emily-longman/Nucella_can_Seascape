# Perform genomic offset analyses with Bio-Oracle data
# Note: prior to running the R script, need to load R and GDAL module on the VACC
# module load R/4.4.1
# module load gdal

# Clear memory
rm(list=ls()) 

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
install.packages(c('data.table', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis', 'terra', 'raster', 'SeqArray', 'poolfstat'))
library(data.table)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(terra)
library(raster)
require(poolfstat)

# Load gradientForest
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
require(gradientForest)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

out_dir <- paste("data/processed/genomic_offset")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

out_fig_dir <- paste("output/figures/genomic_offset")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load SNPs of interest (baypass POD outlier SNPs - 3,095 SNPs SNPs)
baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.01", header=T)

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Metadata
meta <- read.csv("data/raw/pooldata/Populations_metadata.csv", header=T)

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# ================================================================================== #

# Format data

# Create SNP_id column for outlier SNP list
baypass_POD_sig_SNPs <- baypass_POD_sig_SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# For environmental data need a Site-by-Enviro matrix

# Rename "location" column as "sampleId"
names(bio_oracle_sites_2010)[names(bio_oracle_sites_2010) == "location"] <- "Site"

# Rename rows
row.names(bio_oracle_sites_2010) <- bio_oracle_sites_2010$Site

# Subset bio-oracle data so just environmental variables
bio_oracle_sites_2010_sub <- bio_oracle_sites_2010[, c(4:12)]

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Calculate allele frequencies

# Subset pooldata for Baypass outlier SNPs

# Extract SNP info for all SNPs
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  snp.info

# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Filter pooldata for Baypass SNPs of interest
selected_SNPs <- snp.info %>% filter(snp.info$SNP_id %in% baypass_POD_sig_SNPs$SNP_id)
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

# Extract and manipulate coverage for SNPs of interest
coverage %>% as.data.frame -> cov

# Rename columns (19 sites)
names(cov) = c(pooldata_subset@poolnames)

# Transpose coverage
cov_t <- t(cov)

# Extract and manipulate coverage for SNPs of interest
#coverage %>%
#  as.data.frame %>%
#  mutate(SNP_id = snp.info.subset$SNP_id) ->
#  covs.id

# Rename columns (19 sites plus SNP_id)
#names(covs.id) = c(pooldata_subset@poolnames, "SNP_id")

# Restructure data so long format
#reshape2::melt(covs.id, id = "SNP_id", variable.name = "Site", value.name = "COV") -> covs.id.melt

# ================================================================================== #

# Calculate allele frequencies for significant SNPs

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage

# Extract and manipulate allele freq for SNPs of interest
allele_freqs %>% as.data.frame -> afs

# Rename columns (19 sites plus SNP_id)
names(afs) = c(pooldata_subset@poolnames)

# Transpose allele freqs
afs_t <- t(afs)

colnames(afs_t) <- snp.info.subset$SNP_id

# Extract allele freq for SNPs of interest
#allele_freqs %>%
#  as.data.frame %>%
#  mutate(SNP_id = snp.info.subset$SNP_id) ->
#  afs.id

# Rename columns (19 sites and snp_id)
#names(afs.id) = c(pooldata_subset@poolnames, "SNP_id")

# Restructure data so long format
#reshape2::melt(afs.id, id = "SNP_id", variable.name = "Site", value.name = "AF") -> afs.id.melt

# ================================================================================== #

# Sample size of each pool
nSnail=20

# Calcualte mean effective coverage 
nEff <- round((cov_t*2*nSnail)/(cov_t+2*nSnail-1))
# Calculate the effective allele freq
af_nEff <- round((afs_t*nEff)/nEff)

# Check for NAs -- gradient forest can't run if the data includes NAs
which(is.na(af_nEff))

# Join datasets
#left_join(covs.id.melt, afs.id.melt, by = join_by(SNP_id, Site)) -> afs.cov.id

# Calculate mean effective coverage ('nEff')
#afs.cov.id %>% mutate(nEff:=round((COV*2*nSnail)/(COV+2*nSnail-1))) %>%
#  mutate(af_nEff:=round(AF*nEff)/nEff) ->
#  afs.cov.id

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Run gradient forest

nSites <- dim(bio_oracle_sites_2010_sub)[1]
nSpec <- dim(af_nEff)[2]

# Calc maximum number of splits for conditional permutation
lev <- floor(log2(nSites * 0.368/2))
# Note: For conditional permutation, the predictor to be assessed is permuted only within blocks of the dataset defined by splits in the given tree
# on any other predictors correlated above a certain threshold and up to a maximum number of splits set by the maxLevel option 

# Run gradient forest
gf <- gradientForest(cbind(af_nEff, bio_oracle_sites_2010_sub), predictor.vars=colnames(bio_oracle_sites_2010_sub),
                                  response.vars=colnames(af_nEff), ntree=5000, 
                                  maxLevel=lev, trace=T, corr.threshold=0.5)
# Note: Filter by correlation threshold of 0.5
#### Q! Got warnings saying:
# In randomForest.default(x = X, y = spec_vec, maxLevel = maxLevel,  ... :
# The response has five or fewer unique values.  Are you sure you want to do regression?

gf
# Important variables:
# [1] thetao_min   thetao_range thetao_mean  thetao_max   ph_mean

# ================================================================================== #

# Graphing

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
cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(3.1, 1.5, 0.1, 1)))
dev.off()
# Error:
# Error in integrate(approxfun(d, rule = 2), lower = min(d$x), upper = max(d$x)) : 
# roundoff error was detected
# In density.default(splits, weight = w/sum(w), from = rX[1], to = rX[2]) :
# Selecting bandwidth *not* using 'weights'

# Species cumulative plot (For each species shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations)
pdf("output/figures/genomic_offset/species_cumulative_plot.pdf", width = 8, height = 8)
plot(gf, plot.type = "C", imp.vars = most_important, show.overall = F, legend = T, leg.posn = "topleft",
leg.nspecies = 5, cex.lab = 0.7, cex.legend = 0.4, cex.axis = 0.6, line.ylab = 0.9, 
par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))
dev.off()

# Predictor cumulative plot (common.scale=T ensures that plots for all predictors have the same y-scale)
# (For each predictor shows cumulative importance distributions of splits improvement scaled by R2 weighted importance, and standardised by density of observations, averaged over all species.)
pdf("output/figures/genomic_offset/predictor_cumulative_plot.pdf", width = 8, height = 8)
plot(gf, plot.type = "C", imp.vars = most_important, show.species = F, common.scale = T, 
cex.axis = 0.6, cex.lab = 0.7, line.ylab = 0.9, par.args = list(mgp = c(1.5, 0.5, 0), mar = c(2.5, 1, 0.1, 0.5), omi = c(0, 0.3, 0, 0)))
dev.off()

# Fit of Random Forest
pdf("output/figures/genomic_offset/grad_forest_fit.pdf", width = 8, height = 8)
plot(gf, plot.type = "P", show.names = T, horizontal = F, cex.axis = 1, cex.labels = 0.7, line = 2.5)
dev.off()

# ================================================================================== #

# Extract gradient forest results and graph

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
theme(panel.grid.major = element_blank(), legend.position=c(0.8,0.3),legend.key.size = unit(1, 'lines'),
        legend.title = element_blank(), panel.border = element_blank(),axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12), axis.line.x = element_line(colour="black"), 
        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),axis.text.y=element_text(size = 12),axis.title=element_text(size=14,color="black"),
        panel.background = element_blank())+ylab("Accuracy Importance")+ xlab("")
dev.off()

# ================================================================================== #

# Gradient Forest Predictions

# Transform present data
predOUT_present <- predict(gf, bio_oracle_sites_2010_sub)
predOUT_present

# ================================================================================== #

###read in future data
bio_oracle_sites_ssp585_2090 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_ssp585_2090.csv")
rownames(bio_oracle_ssp585_sites) <- rownames(indmeta)




# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Alternate way to get allele frequencies

# Get SNP data

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Filter for sig SNPs

# Filter snp.dt for only significant SNPs
sig_SNPs <- snp.dt %>% filter(snp.dt$SNP_id %in% baypass_POD_sig_SNPs$SNP_id)

# Reset filter
seqResetFilter(genofile)

# Set filter in GDS to only those SNPs
seqSetFilter(genofile, variant.id = sig_SNPs$variant.id)

# Extract sample id and variant id
sample.id = seqGetData(genofile, "sample.id")
variant.id = seqGetData(genofile, "variant.id")

# ================================================================================== #

# Extract allele depth ('ad') of alternate allele for sig SNPs

# Extract allele depth
ad <- seqGetData(genofile, "annotation/format/AD") %>% .$data
# Indexes for alternate allele
even_indexes <- seq(2,ncol(ad),2) # Question: why is there 6193 (ncol(ad)) rather than 6190 (i.e., 3095*2)
# Extract only alternate allele depth
ad_alt <- ad %>% .[,even_indexes]
row.names(ad_alt) = sample.id
#colnames(ad_alt) = variant.id #Doesn't work bc numbers don't align

# Extract total depth ('dp') for sig SNP 
dp <- seqGetData(genofile, "annotation/format/DP")
row.names(dp) = sample.id
colnames(dp) = variant.id

# Transpose and add variant id column
dp %>% t() %>% as.data.frame() %>% mutate(variant.id = variant.id) -> dp.test
# Restructure data so long format
reshape2::melt(dp, id = "variant.id", variable.name = "Site", value.name = "dp") -> dp.melt


# Can't calculate af because ad and dp are not same dimensions
# Calculate allele freq
af <- data.table(ad_alt/dp)

# Merge allele freq table af_i and snp.dt
af_snp <- merge(af, snp.dt, by="variant.id")

# ================================================================================== #




