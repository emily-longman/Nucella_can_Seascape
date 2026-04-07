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
bf.drilling.mean <- foreach(i=0:2, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/phenotypic_data/Baypass_GWAS/NC_pheno_run", i, "_summary_betai_reg.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}
# Change column names
setnames(bf.drilling.mean, "BF(dB)", "bf_db")
# Average BF across replicate runs
bf.drilling.mean.sum <- bf.drilling.mean %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(bf_db_pheno.mean = mean(bf_db), bf_db_pheno.median = median(bf_db),
            bf_db_pheno.var=var(bf_db), eBPis_pheno.mean=mean(eBPis), eBPis_pheno.median=median(eBPis), eBPis_pheno.var=var(eBPis))

# Save
save(bf.drilling.mean.sum, file="data/processed/phenotypic_data/bf.drilling.mean.sum.Rdata")

# Create list of file names
file_names = as.list(dir(path = 'data/processed/phenotypic_data/Baypass_GWAS_POD/', pattern = "*summary_betai_reg.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/phenotypic_data/Baypass_GWAS_POD//', x))))
# Read all the files and add a column with the chunk
bf.POD <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/phenotypic_data/Baypass_GWAS_POD/NC_pheno_POD_run*"))
    # Remove end of chunk name
    tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_betai_reg.out"))
    #Return
    return(tmp)
}
# Change column names
setnames(bf.POD, "BF(dB)", "bf_db")
# Calculate quantiles for each POD
bf.POD.sum <- bf.POD %>% group_by(run) %>% reframe(bf_db = quantile(bf_db, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()
# Average quantiles across POD runs
bf.POD.thr.pheno <- bf.POD.sum %>% group_by(thr) %>% summarize(bf_db.mean=mean(bf_db))

####

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
ggplot(bf.drilling.mean.sum[which(bf.drilling.mean.sum$bf_db_pheno.mean>0),], aes(y=bf_db_pheno.mean, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 20)
dev.off()

# ================================================================================== #
# ================================================================================== #

# Join pheno and pH
bf.ph.mean.sum.pheno <- left_join(bf.ph.mean.sum, bf.drilling.mean.sum, by="MRK")

# Graph phenotypic data vs ph
pdf("output/figures/phenotypic_data/BF_pheno_vs_ph.pdf", width = 8, height = 8)
ggplot(bf.ph.mean.sum.pheno, aes(x=bf_db_pheno.mean, y=bf_db.mean)) +
  labs(x = "BF Pheno", y = "BF pH") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

pdf("output/figures/phenotypic_data/BF_pheno_vs_ph_posBF.pdf", width = 8, height = 8)
ggplot(bf.ph.mean.sum.pheno[which(bf.ph.mean.sum.pheno$bf_db_pheno.mean>0 & bf.ph.mean.sum.pheno$bf_db.mean>0),], 
    aes(x=bf_db_pheno.mean, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF pH") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

# Identify which SNPs beat both thresholds (1 SNPs)
bf.ph.drilling.outliers <- bf.ph.mean.sum.pheno[which(
      bf.ph.mean.sum.pheno$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)] & 
      bf.ph.mean.sum.pheno$bf_db_pheno.mean > bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)]),]
# Write output
write.csv(bf.ph.drilling.outliers, "data/processed/phenotypic_data/bf.ph.drilling.outliers.csv", row.names = F, quote = F)

# ================================================================================== #

# Join pheno and M. tross
bf.Mtross.mean.sum.pheno <- left_join(bf.Mtross.mean.sum, bf.drilling.mean.sum, by="MRK")

# Graph phenotypic data vs M. tross
pdf("output/figures/phenotypic_data/BF_pheno_vs_Mtross.pdf", width = 8, height = 8)
ggplot(bf.Mtross.mean.sum.pheno, aes(x=bf_db_pheno.mean, y=bf_db.mean)) +
  labs(x = "BF Pheno", y = "BF Mtross") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.Mtross$bf_db.mean[which(bf.POD.thr.Mtross$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

pdf("output/figures/phenotypic_data/BF_pheno_vs_Mtross_posBF.pdf", width = 8, height = 8)
ggplot(bf.Mtross.mean.sum.pheno[which(bf.Mtross.mean.sum.pheno$bf_db_pheno.mean>0 & bf.Mtross.mean.sum.pheno$bf_db.mean>0),], 
    aes(x=bf_db_pheno.mean, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF Mtross") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.Mtross$bf_db.mean[which(bf.POD.thr.Mtross$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

# Identify which SNPs beat both thresholds (6 SNPs)
bf.Mtross.drilling.outliers <- bf.Mtross.mean.sum.pheno[which(
      bf.Mtross.mean.sum.pheno$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)] & 
      bf.Mtross.mean.sum.pheno$bf_db_pheno.mean > bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)]),]

# Write output
write.csv(bf.Mtross.drilling.outliers, "data/processed/phenotypic_data/bf.Mtross.drilling.outliers.csv", row.names = F, quote = F)

# ================================================================================== #

# Join pheno and Mcali Thk
bf.McaliIntThk.mean.sum.pheno <- left_join(bf.McaliIntThk.mean.sum, bf.drilling.mean.sum, by="MRK")

# Graph phenotypic data vs M.cali thk
pdf("output/figures/phenotypic_data/BF_pheno_vs_Mcalithk.pdf", width = 8, height = 8)
ggplot(bf.McaliIntThk.mean.sum.pheno, aes(x=bf_db_pheno.mean, y=bf_db.mean)) +
  labs(x = "BF Pheno", y = "BF Mcali_thk") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

pdf("output/figures/phenotypic_data/BF_pheno_vs_Mcalithk_posBF.pdf", width = 8, height = 8)
ggplot(bf.McaliIntThk.mean.sum.pheno[which(bf.McaliIntThk.mean.sum.pheno$bf_db_pheno.mean>0 & bf.McaliIntThk.mean.sum.pheno$bf_db.mean>0),], 
    aes(x=bf_db_pheno.mean, y=bf_db.mean)) + 
  labs(x = "BF Pheno", y = "BF Mcali_thk") +
  geom_point(alpha=0.6) +
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  geom_vline(xintercept=bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)], col="red") +
  theme_classic(base_size = 26)
dev.off()

# Identify which SNPs beat both thresholds (11 SNPs)
bf.McaliIntThk.drilling.outliers <- bf.McaliIntThk.mean.sum.pheno[which(
      bf.McaliIntThk.mean.sum.pheno$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)] & 
      bf.McaliIntThk.mean.sum.pheno$bf_db_pheno.mean > bf.POD.thr.pheno$bf_db.mean[which(bf.POD.thr.pheno$thr==0.999)]),]

# Write output
write.csv(bf.McaliIntThk.drilling.outliers, "data/processed/phenotypic_data/bf.McaliIntThk.drilling.outliers.csv", row.names = F, quote = F)
