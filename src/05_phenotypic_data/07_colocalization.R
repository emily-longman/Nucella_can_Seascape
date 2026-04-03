# Co-localization of BF of GWAS and GEA

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
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)

# ================================================================================== #

# Load data

# Phenotypic GWAS
# Read in SNP data
snp.meta <- read.table("data/processed/phenotypic_data/baypass_input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
# Load phenotypic data
bf.drilling <- read.table("data/processed/phenotypic_data/Baypass_GWAS/NC_pheno_summary_betai_reg.out", header=T)
# Add snpdet
bf.drilling <- cbind(snp.meta, bf.drilling)
# Change column names
setnames(bf.drilling, "BF.dB.", "bf_db")

# Load GEA data
# pH Data
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
# M. tross abundance
load("data/processed/baypass/biotic/bf.Mtross.mean.sum.Rdata")
# M. cali thk
load("data/processed/baypass/biotic/bf.Mcali.IntThk.sum.Rdata")

# Load thresholds
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")
bf.POD.thr.ph <- bf.POD.thr
load("data/processed/baypass/biotic/Mtross_mean_POD_thr.Rdata")
bf.POD.thr.Mtross <- bf.POD.thr
load("data/processed/baypass/biotic/Mcali_IntegratedThk_POD_thr.Rdata")
bf.POD.thr.McaliThk <- bf.POD.thr

# ================================================================================== #

# Graph pheno - only BF > 0
pdf("output/figures/phenotypic_data/BF_pheno_posBF.pdf", width = 15, height = 8)
ggplot(bf.drilling[which(bf.drilling$bf_db>0),], aes(y=bf_db, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.6) + 
  theme_classic(base_size = 20)
dev.off()

# ================================================================================== #
# ================================================================================== #

# Join pheno and pH
bf.ph.mean.sum.pheno <- left_join(bf.ph.mean.sum, bf.drilling, by="MRK")

# Graph phenotypic data vs ph
pdf("output/figures/phenotypic_data/BF_pheno_vs_ph.pdf", width = 8, height = 8)
ggplot(bf.ph.mean.sum.pheno, aes(x=bf_db, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF pH") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  theme_classic(base_size = 20)
dev.off()

# ================================================================================== #

# Join pheno and M. tross
bf.Mtross.mean.sum.pheno <- left_join(bf.Mtross.mean.sum, bf.drilling, by="MRK")

# Graph phenotypic data vs ph
pdf("output/figures/phenotypic_data/BF_pheno_vs_Mtross.pdf", width = 8, height = 8)
ggplot(bf.Mtross.mean.sum.pheno, aes(x=bf_db, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF Mtross") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr.Mtross$bf_db.mean[which(bf.POD.thr.Mtross$thr==0.999)], col="red") +
  theme_classic(base_size = 20)
dev.off()

# ================================================================================== #

# Join pheno and Mcali Thk
bf.McaliIntThk.mean.sum.pheno <- left_join(bf.McaliIntThk.mean.sum, bf.drilling, by="MRK")

# Graph phenotypic data vs ph
pdf("output/figures/phenotypic_data/BF_pheno_vs_Mcalithk.pdf", width = 8, height = 8)
ggplot(bf.McaliIntThk.mean.sum.pheno, aes(x=bf_db, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF Mcali_thk") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  theme_classic(base_size = 20)
dev.off()
