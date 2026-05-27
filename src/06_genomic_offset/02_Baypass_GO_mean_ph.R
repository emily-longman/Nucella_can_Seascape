# Geometric genetic offset

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
    #Return
    return(tmp)
}

# Average across 5 runs
regfile <- regfile_all %>% group_by(COVARIABLE, MRK) %>% summarise(across(everything(), mean))
regfile <-  data.table(regfile)

# Save output
write.table(regfile, "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", row.names=FALSE, col.names=TRUE)
regfile <- read.table("data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", header=T)

# ================================================================================== #

# Compute offset
Ncan_GO <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
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

# Graph GO values per pop

# Color pallet
viridiscolors <- viridis(n=19)

# Order
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

# Redo analyses with scale covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.scaled <- read.table("guide_files/Baypass_ph_mean_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

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
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map.pdf", width = 8, height = 8)s
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 32) + 
  theme(legend.title = element_text(size = 28), legend.text = element_text(size = 20), legend.position = c(0.91, 0.51))
dev.off()
# Graph gGO - alt coloring
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map_alt.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=rev(brewer.pal(6, "BrBG")), name="gGO") +
             coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 32) + 
  theme(legend.title = element_text(size = 28), legend.text = element_text(size = 16), legend.position = c(0.9, 0.51))
dev.off()

