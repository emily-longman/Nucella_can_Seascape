# Geometric genetic offset

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2', 'poolfstat', 'ggplot2', 'RColorBrewer', 'viridis', 'colorspace'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)
library(poolfstat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(colorspace)

# Baypass functions
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")

# ================================================================================== #

# Color palette
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(19))
viridiscolors <- viridis(n=19)

# ================================================================================== #

# Load data

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
# ================================================================================== #

# Analyses with scaled covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.ssp585.scaled <- read.table("guide_files/Baypass_ph_mean_ssp585_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled_ssp585 <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.ssp585.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.ssp585.matrix <- Ncan_GO_scaled_ssp585$go
rownames(go.scaled.ssp585.matrix) <- pooldata@poolnames
colnames(go.scaled.ssp585.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.ssp585.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.ssp585.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.ssp585.output <- left_join(meta, go.scaled.ssp585.output, by="Site")

# Order
go.scaled.ssp585.output$Site <- factor(go.scaled.ssp585.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.ssp585.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
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
pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_GO_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.ssp585.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #

# Extract matrix of RONA estimates between all reference (rows) and target environments (columns)
rona.scaled.ssp585.matrix <- Ncan_GO_scaled_ssp585$rona
rownames(rona.scaled.ssp585.matrix) <- pooldata@poolnames
colnames(rona.scaled.ssp585.matrix) <- pooldata@poolnames

# Extract diagonal
rona.scaled <- diag(rona.scaled.ssp585.matrix)

# Make Site a column
Site <- names(rona.scaled)

# Make table
rona.scaled.ssp585.output <- data.table(Site, rona.scaled)

# Join with metadata
rona.scaled.ssp585.output <- left_join(meta, rona.scaled.ssp585.output, by="Site")

# Order
rona.scaled.ssp585.output$Site <- factor(rona.scaled.ssp585.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_rona.pdf", width = 8, height = 14)
ggplot(rona.scaled.ssp585.output, aes(x = Site, y = rona.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(mycolors)) + ylab("RONA") +
coord_flip() +
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph RONA 
pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_rona_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp585.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_rona_map_grey_shorter.pdf", width = 6.8, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp585.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -112.9), breaks = seq(-125, -112.9, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.745, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

pdf("output/figures/genomic_offset/Baypass_scaled_ssp585_rona_map_grey_wider.pdf", width = 7.05, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp585.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -109.75), breaks = seq(-125, -109.75, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.79, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()


# ================================================================================== #
# ================================================================================== #

# Analyses with scaled covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.ssp370.scaled <- read.table("guide_files/Baypass_ph_mean_ssp370_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled_ssp370 <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.ssp370.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.ssp370.matrix <- Ncan_GO_scaled_ssp370$go
rownames(go.scaled.ssp370.matrix) <- pooldata@poolnames
colnames(go.scaled.ssp370.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.ssp370.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.ssp370.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.ssp370.output <- left_join(meta, go.scaled.ssp370.output, by="Site")

# Order
go.scaled.ssp370.output$Site <- factor(go.scaled.ssp370.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_ssp370_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.ssp370.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
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
pdf("output/figures/genomic_offset/Baypass_scaled_ssp370_GO_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.ssp370.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #

# Extract matrix of RONA estimates between all reference (rows) and target environments (columns)
rona.scaled.ssp370.matrix <- Ncan_GO_scaled_ssp370$rona
rownames(rona.scaled.ssp370.matrix) <- pooldata@poolnames
colnames(rona.scaled.ssp370.matrix) <- pooldata@poolnames

# Extract diagonal
rona.scaled <- diag(rona.scaled.ssp370.matrix)

# Make Site a column
Site <- names(rona.scaled)

# Make table
rona.scaled.ssp370.output <- data.table(Site, rona.scaled)

# Join with metadata
rona.scaled.ssp370.output <- left_join(meta, rona.scaled.ssp370.output, by="Site")

# Order
rona.scaled.ssp370.output$Site <- factor(rona.scaled.ssp370.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_scaled_ssp370_rona.pdf", width = 8, height = 14)
ggplot(rona.scaled.ssp370.output, aes(x = Site, y = rona.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(mycolors)) + ylab("RONA") +
coord_flip() +
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph RONA 
pdf("output/figures/genomic_offset/Baypass_scaled_ssp370_rona_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp370.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA", breaks = c(0.140, 0.145, 0.150)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Analyses with scaled covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.ssp245.scaled <- read.table("guide_files/Baypass_ph_mean_ssp245_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled_ssp245 <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.ssp245.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.ssp245.matrix <- Ncan_GO_scaled_ssp245$go
rownames(go.scaled.ssp245.matrix) <- pooldata@poolnames
colnames(go.scaled.ssp245.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.ssp245.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.ssp245.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.ssp245.output <- left_join(meta, go.scaled.ssp245.output, by="Site")

# Order
go.scaled.ssp245.output$Site <- factor(go.scaled.ssp245.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.ssp245.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
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
pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_GO_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.ssp245.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #

# Extract matrix of RONA estimates between all reference (rows) and target environments (columns)
rona.scaled.ssp245.matrix <- Ncan_GO_scaled_ssp245$rona
rownames(rona.scaled.ssp245.matrix) <- pooldata@poolnames
colnames(rona.scaled.ssp245.matrix) <- pooldata@poolnames

# Extract diagonal
rona.scaled <- diag(rona.scaled.ssp245.matrix)

# Make Site a column
Site <- names(rona.scaled)

# Make table
rona.scaled.ssp245.output <- data.table(Site, rona.scaled)

# Join with metadata
rona.scaled.ssp245.output <- left_join(meta, rona.scaled.ssp245.output, by="Site")

# Order
rona.scaled.ssp245.output$Site <- factor(rona.scaled.ssp245.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_rona.pdf", width = 8, height = 14)
ggplot(rona.scaled.ssp245.output, aes(x = Site, y = rona.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(mycolors)) + ylab("RONA") +
coord_flip() +
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph RONA 
pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_rona_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp245.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA", breaks = c(0.082, 0.084, 0.086)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_rona_map_grey_shorter.pdf", width = 6.8, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp245.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA", breaks = c(0.082, 0.084, 0.086)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -112.9), breaks = seq(-125, -112.9, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.745, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

pdf("output/figures/genomic_offset/Baypass_scaled_ssp245_rona_map_grey_wider.pdf", width = 7.05, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp245.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA", breaks = c(0.082, 0.084, 0.086)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -109.75), breaks = seq(-125, -109.75, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.79, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Analyses with scaled covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.ssp126.scaled <- read.table("guide_files/Baypass_ph_mean_ssp126_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled_ssp126 <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.ssp126.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.ssp126.matrix <- Ncan_GO_scaled_ssp126$go
rownames(go.scaled.ssp126.matrix) <- pooldata@poolnames
colnames(go.scaled.ssp126.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.ssp126.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.ssp126.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.ssp126.output <- left_join(meta, go.scaled.ssp126.output, by="Site")

# Order
go.scaled.ssp126.output$Site <- factor(go.scaled.ssp126.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_ssp126_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.ssp126.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
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
pdf("output/figures/genomic_offset/Baypass_scaled_ssp126_GO_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.ssp126.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #

# Extract matrix of RONA estimates between all reference (rows) and target environments (columns)
rona.scaled.ssp126.matrix <- Ncan_GO_scaled_ssp126$rona
rownames(rona.scaled.ssp126.matrix) <- pooldata@poolnames
colnames(rona.scaled.ssp126.matrix) <- pooldata@poolnames

# Extract diagonal
rona.scaled <- diag(rona.scaled.ssp126.matrix)

# Make Site a column
Site <- names(rona.scaled)

# Make table
rona.scaled.ssp126.output <- data.table(Site, rona.scaled)

# Join with metadata
rona.scaled.ssp126.output <- left_join(meta, rona.scaled.ssp126.output, by="Site")

# Order
rona.scaled.ssp126.output$Site <- factor(rona.scaled.ssp126.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_scaled_ssp126_rona.pdf", width = 8, height = 14)
ggplot(rona.scaled.ssp126.output, aes(x = Site, y = rona.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(mycolors)) + ylab("RONA") +
coord_flip() +
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph RONA 
pdf("output/figures/genomic_offset/Baypass_scaled_ssp126_rona_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp126.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA", breaks = c(0.0282, 0.029, 0.0298)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.75, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Analyses with scaled covariable and original regression coefficients

# Read current and future scaled ph cov file
ph.cov.file.scaled <- read.table("guide_files/Baypass_ph_mean_scaled.txt")
ph.cov.file.future.ssp119.scaled <- read.table("guide_files/Baypass_ph_mean_ssp119_future_scaled.txt")

# Calc GO for scaled pH covar
Ncan_GO_scaled_ssp119 <- compute_genetic_offset(
        beta.coef = NULL, 
        regfile = "data/processed/baypass/abiotic/NC_abiotic_ph_mean_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.ssp119.scaled, scalecov = TRUE, compute.rona = TRUE)

# ================================================================================== #

# Extract matrix of gGO estimates between all reference (rows) and target environments (columns)
go.scaled.ssp119.matrix <- Ncan_GO_scaled_ssp119$go
rownames(go.scaled.ssp119.matrix) <- pooldata@poolnames
colnames(go.scaled.ssp119.matrix) <- pooldata@poolnames

# Extract diagonal
GO.scaled <- diag(go.scaled.ssp119.matrix)

# Make Site a column
Site <- names(GO.scaled)

# Make table
go.scaled.ssp119.output <- data.table(Site, GO.scaled)

# Join with metadata
go.scaled.ssp119.output <- left_join(meta, go.scaled.ssp119.output, by="Site")

# Order
go.scaled.ssp119.output$Site <- factor(go.scaled.ssp119.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_GO.pdf", width = 8, height = 14)
ggplot(go.scaled.ssp119.output, aes(x = Site, y = GO.scaled, fill = Site)) + geom_col() + 
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
pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_GO_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.ssp119.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #

# Extract matrix of RONA estimates between all reference (rows) and target environments (columns)
rona.scaled.ssp119.matrix <- Ncan_GO_scaled_ssp119$rona
rownames(rona.scaled.ssp119.matrix) <- pooldata@poolnames
colnames(rona.scaled.ssp119.matrix) <- pooldata@poolnames

# Extract diagonal
rona.scaled <- diag(rona.scaled.ssp119.matrix)

# Make Site a column
Site <- names(rona.scaled)

# Make table
rona.scaled.ssp119.output <- data.table(Site, rona.scaled)

# Join with metadata
rona.scaled.ssp119.output <- left_join(meta, rona.scaled.ssp119.output, by="Site")

# Order
rona.scaled.ssp119.output$Site <- factor(rona.scaled.ssp119.output$Site, levels=c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL", "PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC"))

# Graph RONA
pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_rona.pdf", width = 8, height = 14)
ggplot(rona.scaled.ssp119.output, aes(x = Site, y = rona.scaled, fill = Site)) + geom_col() + 
scale_fill_manual(values = rev(mycolors)) + ylab("RONA") +
coord_flip() +
theme_bw(base_size = 24) + theme(legend.position="none")
dev.off()

# Graph as map
# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Graph RONA 
pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_rona_map_grey.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp119.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 8) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.76, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_rona_map_grey_shorter.pdf", width = 6.8, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp119.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -112.9), breaks = seq(-125, -112.9, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.745, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

pdf("output/figures/genomic_offset/Baypass_scaled_ssp119_rona_map_grey_wider.pdf", width = 7.05, height = 6.5)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = rona.scaled.ssp119.output, aes(x = Long, y = Lat, fill = rona.scaled), shape = 21, size = 6) + 
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="RONA") +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -109.75), breaks = seq(-125, -109.75, by = 5)) + ylim(32.5, 46.25) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.79, 0.77), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Put all in one dataframe to graph with similar scale

# Add scenario
rona.scaled.ssp585.output$scenario <- "SSP5-8.5"
rona.scaled.ssp370.output$scenario <- "SSP3-7.0"
rona.scaled.ssp245.output$scenario <- "SSP2-4.5"
rona.scaled.ssp126.output$scenario <- "SSP1-2.6"
rona.scaled.ssp119.output$scenario <- "SSP1-1.9"

# Join
rona.all <- bind_rows(rona.scaled.ssp119.output, rona.scaled.ssp126.output, rona.scaled.ssp245.output, rona.scaled.ssp370.output, rona.scaled.ssp585.output)

# Set colors
ssp_colors = sequential_hcl(6, palette = "Blues 3")

# Graph RONA for all SSP scenarios
pdf("output/figures/genomic_offset/Baypass_scaled_sspALL_rona.pdf", width = 5.75, height = 6)
ggplot(rona.all, aes(x = Lat, y = rona.scaled, fill = scenario))+
  geom_point(shape = 21, size = 8) + 
  scale_fill_manual(values = rev(ssp_colors[2:6])) +
  scale_x_continuous(limits = c(34, 45.5), breaks = c(36, 40, 44)) +
  scale_y_continuous(limits = c(-0.01, 0.21), breaks = c(0, 0.1, 0.2)) +
  labs(x = "Latitude", y = "RONA") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()

# Graph RONA for all SSP scenarios - lines
pdf("output/figures/genomic_offset/Baypass_scaled_sspALL_rona_lines.pdf", width = 5.75, height = 6)
ggplot(rona.all, aes(x = Lat, y = rona.scaled, color = scenario))+
  geom_line(linewidth = 4) + 
  scale_color_manual(values = rev(ssp_colors[1:5])) +
  scale_x_continuous(limits = c(34.5, 45.25), breaks = c(36, 40, 44)) +
  scale_y_continuous(limits = c(-0.01, 0.21), breaks = c(0, 0.1, 0.2)) +
  labs(x = "Latitude", y = "RONA") +
  theme_linedraw(base_size = 30) + theme(legend.position = "none")
dev.off()
