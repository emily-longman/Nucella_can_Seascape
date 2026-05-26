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

# ================================================================================== #

# Load regfile and summarise

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/abiotic/ph_mean_scaled/', pattern = "*summary_betai_reg.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/abiotic/ph_mean_scaled/', x))))

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
write.table(regfile, "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_scaled_run_all_summary_betai_reg.out", row.names=FALSE, col.names=TRUE)

# ================================================================================== #
# ================================================================================== #

# Redo analyses and scale covariable - note I already did this and used these values for the run of Baypass that I am using above

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
        regfile = "data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_scaled_run_all_summary_betai_reg.out", 
        covfile = "guide_files/Baypass_ph_mean_scaled.txt",
        newenv = ph.cov.file.future.scaled, scalecov = TRUE, compute.rona = TRUE)

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

# Color pallet
viridiscolors <- viridis(n=19)

# Graph GO
pdf("output/figures/genomic_offset/Baypass_scaled_GO_newrun.pdf", width = 8, height = 14)
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
pdf("output/figures/genomic_offset/Baypass_scaled_GO_map_newrun.pdf", width = 8, height = 8)
ggplot(data = west_coast) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = go.scaled.output, aes(x = Long, y = Lat, fill = GO.scaled), shape = 21, size = 8) + 
  #scale_fill_gradient(low = "cyan1", high = "gray27") + 
  #scale_fill_viridis(option="viridis", direction = -1) +
  scale_fill_gradientn(colours=brewer.pal(9, "RdPu"), name="gGO", breaks = c(0.085, 0.090, 0.095)) +
  coord_fixed(1.3) +
  xlim(c(-125, -114)) +
  xlab("Longitude") + ylab("Latitude") + theme_classic(base_size = 27) + 
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.9, 0.52))
dev.off()
