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
out_fig_dir <- paste("output/figures/baypass")
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

# Baypass in standard covariate mode

# Read all files and add variable name
baypass.bf <- foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/abiotic/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_betai_reg.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_betai <- read.table(full_path, header=T)

    # Add variable column
    baypass_betai <- baypass_betai %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_betai_pos <- cbind(snp.meta, baypass_betai)
}

# Graph BF of individual variables
pdf("output/figures/baypass/baypass_BF_ph_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "ph_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/baypass/baypass_BF_thetao_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "thetao_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/baypass/baypass_BF_chl_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "chl_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/baypass/baypass_BF_so_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "so_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

# Read in files and graph - abiotic
foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/abiotic/", w, sep = "")
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

    # Graph XtX for var
    pdf(paste0("output/figures/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for var
    pdf(paste0("output/figures/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
    par(mar=c(5,5,4,1)+.1) #Adjust margins
    plot(baypass_xtx$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
    abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
    abline(h=-log10(0.05/dim(baypass_xtx)[1]), lty=1, col="red") # Bonferroni threshold
    dev.off()

}

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Read all files and add variable name
baypass.bf <- foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/abiotic_aux/", w, sep = "")
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
    path <- paste("data/processed/baypass/biotic/Mcali/", w, sep = "")
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
pdf("output/figures/baypass/baypass_BF.pdf", width = 8, height = 8)
ggplot(baypass.bf, aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  facet_wrap(~variable) +
  theme_bw() + theme(legend.position = "none") + theme_classic(base_size = 20)
dev.off()

# Graph BF of individual variables
pdf("output/figures/baypass/baypass_BF_ph_mean_baseR.pdf", width = 8, height = 8)
plot(baypass.bf[which(baypass.bf$variable == "ph_mean"),]$BF.dB., xlab="SNP",ylab="BFaux (in dB)")
dev.off()

pdf("output/figures/baypass/baypass_BF_ph_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "ph_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/baypass/baypass_BF_thetao_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "thetao_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

pdf("output/figures/baypass/baypass_BF_chl_mean.pdf", width = 8, height = 8)
ggplot(baypass.bf[which(baypass.bf$variable == "chl_mean"),], aes(y=BF.dB., x=chr)) + 
  labs(x = "Position", y = "BFaux (in dB)") +
  geom_point(alpha=0.8) + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()

# ================================================================================== #

# Explore patterns in bayes factors

# pH mean
bf.ph.mean <- baypass.bf[which(baypass.bf$variable == "ph_mean"),]
bf.ph.mean.extreme <- bf.ph.mean[which(bf.ph.mean$BF.dB.>50),]
bf.ph.mean.extreme <- bf.ph.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.ph.mean.extreme, "data/processed/baypass/bf.ph.mean.extremes.csv", row.names = F, quote = F)

# Temp mean
bf.thetao.mean <- baypass.bf[which(baypass.bf$variable == "thetao_mean"),]
bf.thetao.mean.extreme <- bf.thetao.mean[which(bf.thetao.mean$BF.dB.>50),]
bf.thetao.mean.extreme <- bf.thetao.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.thetao.mean.extreme, "data/processed/baypass/bf.thetao.mean.extremes.csv", row.names = F, quote = F)

# Temp min
bf.thetao.min <- baypass.bf[which(baypass.bf$variable == "thetao_min"),]
bf.thetao.min.extreme <- bf.thetao.min[which(bf.thetao.min$BF.dB.>50),]
bf.thetao.min.extreme <- bf.thetao.min.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.thetao.min.extreme, "data/processed/baypass/bf.thetao.min.extremes.csv", row.names = F, quote = F)

# Chl mean
bf.chl.mean <- baypass.bf[which(baypass.bf$variable == "chl_mean"),]
bf.chl.mean.extreme <- bf.chl.mean[which(bf.chl.mean$BF.dB.>50),]
bf.chl.mean.extreme <- bf.chl.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.chl.mean.extreme, "data/processed/baypass/bf.chl.mean.extreme.csv", row.names = F, quote = F)

# Note: All have a max Bayes Factor value of 52.96447989

# No overlap between the SNPs in the ph mean and thetao mean nor
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.thetao.mean.extreme$SNP_id)
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.chl.mean.extreme$SNP_id)
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.thetao.min.extreme$SNP_id)

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# ================================================================================== #
# ================================================================================== #

# Read in files and graph - abiotic
foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/abiotic/", w, sep = "")
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

    # Graph XtX for var
    pdf(paste0("output/figures/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for var
    pdf(paste0("output/figures/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
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
    path <- paste("data/processed/baypass/biotic/Mcali/", w, sep = "")
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
    pdf(paste0("output/figures/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for ph_mean
    pdf(paste0("output/figures/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
    par(mar=c(5,5,4,1)+.1) #Adjust margins
    plot(baypass_xtx$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
    abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
    abline(h=-log10(0.05/dim(baypass_xtx)[1]), lty=1, col="red") # Bonferroni threshold
    dev.off()

}






