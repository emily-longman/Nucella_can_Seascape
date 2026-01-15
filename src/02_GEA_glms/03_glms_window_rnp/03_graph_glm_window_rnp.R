# Graph window analysis

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Merge GLM windows

# Create list of file names
path <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_ph_mean/")
file_names = as.list(dir(path = path, pattern = "glm_window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_ph_mean/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
win.out =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(win.out)

# ================================================================================== #

# Load windows
load("data/processed/GEA/glms/glms_window_summary/windows.RData")

# ================================================================================== #

# Summarize model
win.ph.mean <- win.out %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, win, chr, pos_mean) %>%
  summarise(uci = quantile(rnp.binom.p.0.05, 0.05)) %>%
  mutate(model = "ph_mean")

# ================================================================================== #

# Graph rnp p

# Create unique Chromosome number
chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = chr.unique))
win.ph.mean.chr.unique <- unique(win.ph.mean$chr)
win.ph.mean$chr.unique <- as.numeric(factor(win.ph.mean$chr, levels = win.ph.mean.chr.unique))

# Graph rnp geompoint
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_ph_mean_rnp_0.05_geompoint_real.pdf", width = 8, height = 8)
ggplot(win.out[which(win.out$perm==0),], aes(y=-log10(rnp.binom.p.0.05), x=chr.unique)) + 
  geom_point(alpha=0.8, size=1.3) + 
  facet_wrap(~variable) + 
  theme_bw() + theme(legend.position = "none")
dev.off()

# Graph rnp geomline
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_ph_mean_rnp_0.05_geomline_real.pdf", width = 8, height = 8)
ggplot(win.out[which(win.out$perm==0),], aes(y=-log10(rnp.binom.p.0.05), x=chr.unique)) + 
  geom_line( ) + 
  facet_wrap(~variable) +
  theme_bw() + theme(legend.position = "none")
dev.off()

# Graph permutations
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_ph_mean_uci.pdf", width = 14, height = 8)
ggplot(win.ph.mean, aes(y=-log10(uci), x=chr.unique, col=data_type)) + 
  geom_line() + ylim(0,250) +
  scale_color_manual(values = c("#a1c8cf", "black")) +
  theme_bw()
dev.off()

#pdf("output/figures/GEA/glms/glms_window_summary/glm_window_rnp_0.05.pdf", width = 8, height = 8)
#ggplot(win.out, aes(y=-log10(rnp.binom.p.0.05), x=chr.unique, col=factor(perm))) + 
#  geom_point(alpha=0.8, size=1.3) + 
#  scale_color_manual(values = c("black", rep("red", 100))) +
#  facet_wrap(~variable) +
#  theme_bw() + theme(legend.position = "none")
#dev.off()


# ================================================================================== #
# ================================================================================== #

# Merge GLM windows

# Create list of file names
path <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_Mtross_mean/")
file_names = as.list(dir(path = path, pattern = "glm_window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_Mtross_mean/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
win.out =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(win.out)

# ================================================================================== #

# Summarize model
win.Mtross.mean <- win.out %>%
  mutate(data_type = case_when(perm == 0 ~ "real",
                               perm != 0 ~ "perm")) %>%
  group_by(data_type, win, chr, pos_mean) %>%
  summarise(uci = quantile(rnp.binom.p.0.05, 0.05)) %>%
  mutate(model = "Mtross_mean")

  # ================================================================================== #

# Graph rnp p

# Create unique Chromosome number
chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = chr.unique))
win.Mtross.mean.chr.unique <- unique(win.Mtross.mean$chr)
win.Mtross.mean$chr.unique <- as.numeric(factor(win.Mtross.mean$chr, levels = win.Mtross.mean.chr.unique))

# Graph rnp geompoint
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_Mtross_mean_rnp_0.05_geompoint_real.pdf", width = 8, height = 8)
ggplot(win.out[which(win.out$perm==0),], aes(y=-log10(rnp.binom.p.0.05), x=chr.unique)) + 
  geom_point(alpha=0.8, size=1.3) + 
  facet_wrap(~variable) + 
  theme_bw() + theme(legend.position = "none")
dev.off()

# Graph rnp geomline
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_Mtross_mean_rnp_0.05_geomline_real.pdf", width = 8, height = 8)
ggplot(win.out[which(win.out$perm==0),], aes(y=-log10(rnp.binom.p.0.05), x=chr.unique)) + 
  geom_line( ) + 
  facet_wrap(~variable) +
  theme_bw() + theme(legend.position = "none")
dev.off()

# Graph permutations
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_Mtross_mean_uci.pdf", width = 14, height = 8)
ggplot(win.Mtross.mean, aes(y=-log10(uci), x=chr.unique, col=data_type)) + 
  geom_line() + ylim(0,250) +
  scale_color_manual(values = c("#a1c8cf", "black")) +
  theme_bw()
dev.off()