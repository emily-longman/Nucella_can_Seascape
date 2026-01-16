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
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load POD threshold
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")

# ================================================================================== #
 
# Read in Baypass pH mean files

# Load BF output for 5 replicate Baypass runs
bf.ph.mean <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run", i, "_summary_betai_reg.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Change column names
setnames(bf.ph.mean, "BF(dB)", "bf_db")

# ================================================================================== #

# Average BF across replicate runs
bf.ph.mean.sum <- bf.ph.mean %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(bf_db.mean = mean(bf_db), bf_db.median = median(bf_db),
            bf_db.var=var(bf_db), eBPis.mean=mean(eBPis), eBPis.median=median(eBPis), eBPis.var=var(eBPis))

# ================================================================================== #

# Graph BF with POD threshold
pdf("output/figures/baypass/baypass_BF_ph_mean_rep.means.pdf", width = 12, height = 8)
ggplot(bf.ph.mean.sum, aes(y=bf_db.mean, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))
dev.off()

# Graph BF with POD threshold - without BF < 0 
pdf("output/figures/baypass/baypass_BF_ph_mean_5rep_posBF.pdf", width = 12, height = 8)
ggplot(bf.ph.mean.sum[which(bf.ph.mean.sum$bf_db.mean>0),], aes(y=bf_db.mean, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red") +
  theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))
dev.off()

# ================================================================================== #

# Identify patterns in bayes factors

# pH mean - 1,673 SNPs with BF > threshold
bf.ph.mean.sum.outliers <- bf.ph.mean.sum[which(bf.ph.mean.sum$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)]),]
bf.ph.mean.sum.outliers <- bf.ph.mean.sum.outliers %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.ph.mean.sum.outliers, "data/processed/baypass/bf.ph.mean.sum.outliers.csv", row.names = F, quote = F)

