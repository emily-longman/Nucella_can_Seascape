# Use the PODs to calibrate Bayes Factor

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
install.packages(c('data.table', 'foreach', 'tidyverse', 'ggplot2'))
library(data.table)
library(foreach)
library(tidyverse)
library(ggplot2)

# ================================================================================== #

# Generate Folders and files

# Make output directories
data_processed_outlier="data/processed/baypass"
if (!dir.exists(data_processed_outlier)) {dir.create(data_processed_outlier)}
data_processed_outlier="data/processed/baypass/input_files"
if (!dir.exists(data_processed_outlier)) {dir.create(data_processed_outlier)}

# ================================================================================== #
# ================================================================================== #

# Load POD files - mean pH

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/abiotic/ph_mean_POD/', pattern = "*summary_betai_reg.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/abiotic/ph_mean_POD/', x))))

# Read all the files and add a column with the chunk
bf.POD <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/baypass/abiotic/ph_mean_POD/NC_abiotic_ph_mean_POD_run*"))
    # Remove end of chunk name
    tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_betai_reg.out"))
    #Return
    return(tmp)
}

# Change column names
setnames(bf.POD, "BF(dB)", "bf_db")

# ================================================================================== #

# Calculate quantiles for each POD
bf.POD.sum <- bf.POD %>% group_by(run) %>% reframe(bf_db = quantile(bf_db, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()

# Average quantiles across POD runs
bf.POD.thr <- bf.POD.sum %>% group_by(thr) %>% summarize(bf_db.mean=mean(bf_db))

# Graph BF of one POD to see distribution with threshold
pdf("output/figures/baypass/baypass_ph_mean_BF_POD1.pdf", width = 8, height = 8)
ggplot(bf.POD[which(bf.POD$run==1),], aes(y=bf_db, x=MRK)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.7) + 
  theme_classic(base_size = 20) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red") + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

# Save
save(bf.POD.thr, file="data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")

# ================================================================================== #
# ================================================================================== #

# Load POD files - M californianus integrated shell thk

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/biotic/Mcali_IntegratedThk_POD/', pattern = "*summary_betai_reg.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/biotic/Mcali_IntegratedThk_POD/', x))))

# Read all the files and add a column with the chunk
bf.POD <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/baypass/biotic/Mcali_IntegratedThk_POD/NC_biotic_Mcali_IntegratedThk_POD_run*"))
    # Remove end of chunk name
    tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_betai_reg.out"))
    #Return
    return(tmp)
}

# Change column names
setnames(bf.POD, "BF(dB)", "bf_db")

# ================================================================================== #

# Calculate quantiles for each POD
bf.POD.sum <- bf.POD %>% group_by(run) %>% reframe(bf_db = quantile(bf_db, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()

# Average quantiles across POD runs
bf.POD.thr <- bf.POD.sum %>% group_by(thr) %>% summarize(bf_db.mean=mean(bf_db))

# Graph BF of one POD to see distribution with threshold
pdf("output/figures/baypass/baypass_Mcali_IntegratedThk_BF_POD1.pdf", width = 8, height = 8)
ggplot(bf.POD[which(bf.POD$run==1),], aes(y=bf_db, x=MRK)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.7) + 
  theme_classic(base_size = 20) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red") + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

# Save
save(bf.POD.thr, file="data/processed/baypass/biotic/Mcali_IntegratedThk_POD_thr.Rdata")

