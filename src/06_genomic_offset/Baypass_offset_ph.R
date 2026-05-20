# Try genomic offset with baypass

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
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)
library(poolfstat)
library(RColorBrewer)
library(viridis)

# Baypass functions
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")

# ================================================================================== #

# Load data

# Load ph var in future (2090 ssp 585)
ph.cov.file.future <- read.table("guide_files/Baypass_ph_mean_future.txt")

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Load metadata
meta <- read.csv("guide_files/Populations_metadata.csv", header=T)

# Load regression file (test with just one file)
#regfile_run1 <- read.table("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run1_summary_betai_reg.out", header=T)

# ================================================================================== #

# Load regfile and summarise

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/abiotic/ph_mean/', pattern = "*summary_betai_reg.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/abiotic/ph_mean/', x))))

# Read all the files and add a column with the chunk
regfile_all <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    #tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run*"))
    # Remove end of chunk name
    #tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_betai_reg.out"))
    #Return
    return(tmp)
}

# Average across 5 runs
regfile <- regfile_all %>% group_by(COVARIABLE, MRK) %>% summarise(across(everything(), mean))
regfile <-  data.table(regfile)

# Save output
write.table(regfile, "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", row.names=FALSE, col.names=TRUE)

# ================================================================================== #

# Compute offset
Ncan_GO <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean.txt",
        newenv = ph.cov.file.future, scalecov = FALSE, compute.rona = TRUE)

# Notes:
# beta.coef - matrix of reg coef; if NULL, then provide BayPass output file with regfile argument
# regfile - BayPass output file with estimates of the regression coefficients
# scalecov - if TRUE all covariables are scaled with respect to mean and variance of original covariable value
# compute.rona - if true, calculate the RONA statistic from Rellstab et al. 2016
#   Rona is  related to the square root of the geometric GO (times √2/π)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.matrix <- Ncan_GO$go
rownames(go.matrix) <- pooldata@poolnames
colnames(go.matrix) <- pooldata@poolnames

# Extract diagonal
GO <- diag(go.matrix)

# Make Site a column
Site <- names(GO)

# Make table
go.output <- data.table(Site, GO)

# Join with metadata
go.output <- left_join(meta, go.output, by="Site")

# ================================================================================== #

# Color pallet
viridiscolors <- viridis(n=19)

# Order
#go.output$Site <- factor(go.output$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
go.output$Site <- factor(go.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_GO.pdf", width = 8, height = 14)
ggplot(go.output, aes(x = Site, y = GO, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(viridiscolors)) + ylab("gGO") +
coord_flip() + 
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# ================================================================================== #

# Extract rona matrix of GO estimates between all reference (rows) and target environments (columns)
rona.matrix <- Ncan_GO$rona
rownames(rona.matrix) <- pooldata@poolnames
colnames(rona.matrix) <- pooldata@poolnames

# Extract diagonal
RONA <- diag(rona.matrix)

# Make Site a column
Site <- names(RONA)

# Make table
RONA.output <- data.table(Site, RONA)

# Join with metadata
RONA.output <- left_join(meta, RONA.output, by="Site")

# Order
RONA.output$Site <- factor(RONA.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_RONA.pdf", width = 8, height = 14)
ggplot(RONA.output, aes(x = Site, y = RONA, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(viridiscolors)) + 
coord_flip() + 
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Redo analyses and scale covariable

# Read ph cov file
ph.cov.file <- read.table("guide_files/Baypass_ph_mean.txt")

# Calc mean
ph.mean <- mean(t(ph.cov.file))
# Calc sd
ph.sd <- sd(t(ph.cov.file))

# Scale
ph.cov.file.scaled <- (ph.cov.file - ph.mean)/ph.sd
# Write file
write.table(ph.cov.file.scaled, "guide_files/Baypass_ph_mean_scaled.txt", col.names=F, row.names=F)

# Scale future
ph.cov.file.future.scaled <- (ph.cov.file.future - ph.mean)/ph.sd

# Calc GO for scaled pH covar
Ncan_GO_scaled <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future, scalecov = TRUE, compute.rona = TRUE)

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.matrix <- Ncan_GO_scaled$go
rownames(go.scaled.matrix) <- pooldata@poolnames
colnames(go.scaled.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.output <- left_join(meta, go.scaled.output, by="Site")

# Order
go.scaled.output$Site <- factor(go.scaled.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(viridiscolors)) + ylab("gGO scaled") +
coord_flip() + 
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph gGO
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  #scale_fill_viridis(option="viridis", direction = -1) +
  scale_fill_gradientn(colours=rev(brewer.pal(6, "BrBG")), name="gGO") +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 27) + 
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 16), legend.position = c(0.98, 0.52))
dev.off()

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Perform analyses on subset of SNPs - just pos selected for pH - NOT WORKING

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load and merge xtx data
NC.xtx <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/xtx/NC_run", i, "_summary_pi_xtx.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}
colnames(NC.xtx) <- c("chr", "pos", "allele1", "allele2", "MRK", "M_P", "SD_P", "M_XtX", "SD_XtX", "XtXst", "log10.1.pval.", "rep")
# Average xtx across replicate runs
NC.xtx.sum <- NC.xtx %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(XtXst.mean = mean(XtXst), XtXst.median = median(XtXst), log10.1.pval.mean = mean(log10.1.pval.), log10.1.pval.median = median(log10.1.pval.))


# Calculate mean xtx
mean.xtx <- mean(NC.xtx.sum$XtXst.mean)
# Standardize xtx
NC.xtx.sum <- NC.xtx.sum %>% mutate(XtXst.mean.standardize = XtXst.mean - mean.xtx)
# Identify SNPs with XtXst higher than the mean and are significant
NC.xtx.sum.pos <- NC.xtx.sum %>% filter(XtXst.mean.standardize > 0 & log10.1.pval.mean > -log10(0.001))
# Identify indices of SNPs
indices <- which(NC.xtx.sum$XtXst.mean.standardize > 0 & NC.xtx.sum$log10.1.pval.mean > -log10(0.001))

# Load SNPs of interest - SNPs undergoing positive selection (baypass POD outlier SNPs - 8914 SNPs SNPs)
#load("data/processed/genomic_offset/NC.xtx.sum.pos.Rdata")

#####
# Load Baypass pH output
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")
# Identify SNPs which beat POD threshold
bf.ph.mean.sum.outliers <- bf.ph.mean.sum[which(bf.ph.mean.sum$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)]),]

# Identify indices of SNPs
bf.ph.mean.sum.outliers.MRK <- bf.ph.mean.sum$MRK[which(bf.ph.mean.sum$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)])]

# Compute offset
Ncan_GO_pos <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean.txt",
        newenv = ph.cov.file.future, scalecov = FALSE,
        candidate.snp = bf.ph.mean.sum.outliers.MRK)
# Get error: "Error in eigen(BtB, symmetric = TRUE) : non-square matrix in 'eigen'""


# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

install.packages(c("terra", "geodata", "fields", "maps", "LEA"))

library(terra)
library(geodata)
library(fields)
library(maps)
library(LEA)