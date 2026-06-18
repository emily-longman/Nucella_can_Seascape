# Summarize corrected af values for each pop

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)

# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Subset pops and cut ARA
pooldata.subset.18pop <- pooldata.subset(pooldata, pool.index=c(1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
pops <- data.frame(Site = pooldata.subset.18pop@poolnames, POP = 1:length(pooldata.subset.18pop@poolnames))

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
snp.meta$MRK <- c(1:dim(snp.meta)[1])
snp.meta <- snp.meta %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Read in Baypass mean pH  files

# Load BF output for 5 replicate Baypass runs
baypass.Mcali <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/biotic/Mcali_IntegratedThk/NC_biotic_Mcali_IntegratedThk_run", i, "_summary_yij_pij.out", sep=""))
    tmp[,rep:=i]
    return(tmp)
}

# ================================================================================== #

# Average  across replicate runs
# M_P is the mean of the posterior dist of the aij param, which is closely related to the freq of the reference allele
# M_Y is the posterior mean of allele counts)
baypass.Mcali.sum <- baypass.Mcali %>% group_by(POP, MRK) %>% 
    reframe(M_Y.mean = mean(M_Y), M_P.mean = mean(M_P))

# Join pops with Baypass output
baypass.Mcali.sum <- left_join(baypass.Mcali.sum, pops, by = "POP")

# Join with metadata
baypass.Mcali.sum <- left_join(baypass.Mcali.sum, snp.meta, by = "MRK")

# ================================================================================== #

# Save output
save(baypass.Mcali.sum, file = "data/processed/baypass/biotic/baypass.Mcali.sum.RData")
