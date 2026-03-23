# Create windows

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'doMC'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/baypass/window_summary/window_chunk_analysis_Mcali_IntThk/")
file_names = as.list(dir(path = path, pattern = "window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/baypass/window_summary/window_chunk_analysis_Mcali_IntThk/"), x))))

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

# Save merged data
save(win.out, file = "data/processed/baypass/window_summary/window_analysis_Mcali_IntThk.RData")

# ================================================================================== #

# Use the POD threshold to come up with p-val

# Load mean bf data from 5 baypass runs
load("data/processed/baypass/biotic/bf.Mcali.IntThk.sum.Rdata")
# Load POD thresholds
load("data/processed/baypass/biotic/Mcali_IntegratedThk_POD_thr.Rdata")

# Create the ECDF (empirical cumulative distribution functio) function
my_ecdf <- ecdf(bf.McaliIntThk.mean.sum$bf_db.mean)

# Find the probability for a given value
probability <- my_ecdf(bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)])

# Calc pr.i as the opposite of the probability
pr.i <- c(1-probability)

# ================================================================================== #

# Create unique Chromosome number
win.out.chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = win.out.chr.unique))

# Graph rnp p

# Graph rnp geompoint
pdf("output/figures/baypass/window_summary/baypass_window_Mcali_IntThk_rnpPOD_geompoint.pdf", width = 12, height = 6)
ggplot(win.out, aes(y=-log10(rnp.binom.POD), x=chr.unique)) + 
  geom_point(alpha=0.8, size=1.6) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=26) + theme(legend.position = "none")
dev.off()

# Graph rnp geomline
pdf("output/figures/baypass/window_summary/baypass_window_Mcali_IntThk_rnpPOD_geomline.pdf", width = 12, height = 6)
ggplot(win.out, aes(y=-log10(rnp.binom.POD), x=chr.unique)) + 
  geom_line( ) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=26) + theme(legend.position = "none")
dev.off()

# ================================================================================== #

# Extract outliers
win.out.outliers <- win.out %>% filter(-log10(rnp.binom.POD) > -log10(0.01/length(rnp.binom.POD)))

# Save outliers
write.csv(win.out.outliers, "data/processed/baypass/window_summary/window_analysis_Mcali_IntThk_outliers.csv", row.names=FALSE)

