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
out_fig_dir <- paste("output/figures/GEA/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load Data

# Read in SNP data
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Ecological variables
abiotic_variables <- read.table("guide_files/Seascape_vars_names.txt")
abiotic_var <- abiotic_variables$V1[1:9]
Mcali_variables <- read.table("guide_files/Mcali_vars_names.txt")
Mcali_var <- Mcali_variables$V1[1:4]

# ================================================================================== #

# Read all files and add variable name
baypass.bf <- foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/GEA/baypass/abiotic/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_betai.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_betai <- read.table(full_path, header=T)

    # Add variable column
    baypass_betai <- baypass_betai %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_betai_pos <- cbind(snp.meta, baypass_betai)
}

# Read all files and add variable name
baypass.bf <- foreach(w=Mcali_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/GEA/baypass/biotic/Mcali/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_betai.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_betai <- read.table(full_path, header=T)

    # Add variable column
    baypass_betai <- baypass_betai %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_betai_pos <- cbind(snp.meta, baypass_betai)
}

# ================================================================================== #

# Graph Baypass output

# Graph bayes factor - This is too big.....
pdf("output/figures/GEA/baypass/glm_baypass_BF.pdf", width = 8, height = 8)
ggplot(baypass.bf, aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  facet_wrap(~variable) +
  theme_bw() + theme(legend.position = "none") + theme_classic(base_size = 20)
dev.off()

# Graph BF of individual variables
pdf("output/figures/GEA/baypass/glm_baypass_BF_ph_mean.pdf", width = 8, height = 8)
plot(baypass.bf[which(baypass.bf$variable == "ph_mean"),]$BF.dB., xlab="SNP",ylab="BFmc (in dB)")
dev.off()

pdf("output/figures/GEA/baypass/glm_baypass_BF_ph_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "ph_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/GEA/baypass/glm_baypass_BF_chl_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "chl_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

# ================================================================================== #
# ================================================================================== #

# Read in files and graph - abiotic
foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/GEA/baypass/abiotic/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_summary_pi_xtx.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_xtx <- read.table(full_path, header=T)

    # Add variable column
    baypass_xtx <- baypass_xtx %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_xtx <- cbind(snp.meta, baypass_xtx)

    # Graph XtX for ph_mean
    pdf(paste0("output/figures/GEA/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for ph_mean
    pdf(paste0("output/figures/GEA/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
    par(mar=c(5,5,4,1)+.1) #Adjust margins
    plot(baypass_xtx$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
    abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
    abline(h=-log10(0.05/dim(baypass_xtx)[1]), lty=1, col="red") # Bonferroni threshold
    dev.off()

}

# Read in files and graph - biotic
foreach(w=Mcali_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/GEA/baypass/biotic/Mcali/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_summary_pi_xtx.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_xtx <- read.table(full_path, header=T)

    # Add variable column
    baypass_xtx <- baypass_xtx %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_xtx <- cbind(snp.meta, baypass_xtx)

    # Graph XtX for ph_mean
    pdf(paste0("output/figures/GEA/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for ph_mean
    pdf(paste0("output/figures/GEA/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
    par(mar=c(5,5,4,1)+.1) #Adjust margins
    plot(baypass_xtx$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
    abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
    abline(h=-log10(0.05/dim(baypass_xtx)[1]), lty=1, col="red") # Bonferroni threshold
    dev.off()

}






