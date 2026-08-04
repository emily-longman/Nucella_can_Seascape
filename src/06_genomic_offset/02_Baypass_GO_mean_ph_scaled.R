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
library(colorspace)

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

# Save output
save(go.scaled.output, file = "data/processed/genomic_offset/Nucella_gGO.Rdata")

# ================================================================================== #

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
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 9) + 
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32, 48) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 32) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 28), legend.text = element_text(size = 20), legend.position = c(0.818, 0.51))
dev.off()

# Graph gGO - alt coloring
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map_alt2.pdf", width = 8, height = 8.39)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  #scale_fill_gradientn(colours=brewer.pal(6, "YlGn"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_continuous_sequential(palette = "Purples 2", name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.785, 0.84), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# Graph gGO - alt coloring Grey
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map_alt3_shorter.pdf", width = 8, height = 7)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 7) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + ylim(32.5, 46.5) +
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.745, 0.8), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

pdf("output/figures/genomic_offset/Baypass_scaled_GO_map_alt3_taller.pdf", width = 8, height = 11)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 11) + 
  #scale_fill_viridis(option="rocket", breaks = c(0.085, 0.090, 0.095)) +
  scale_fill_gradientn(colours=brewer.pal(9, "Greys"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3), expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(32, 46.5), breaks = seq(32, 44, by = 4), expand = expansion(mult = c(0.01, 0.01))) +
  #scale_x_continuous(limits = c(-125, -114.1), breaks = seq(-125, -114.1, by = 3)) + 
  #xlim(c(-125, -112.5)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + 
  theme(
    panel.grid.major = element_blank(), # Removes major grid lines
    panel.grid.minor = element_blank(), # Removes minor grid lines
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) + # Keeps outer rectangle
  theme(legend.title = element_text(size = 26), legend.text = element_text(size = 20), legend.position = c(0.83, 0.87), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

