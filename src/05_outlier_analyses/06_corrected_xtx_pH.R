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

pops <- data.frame(Site = pooldata@poolnames, POP = 1:length(pooldata@poolnames))

# ================================================================================== #

# Read in Baypass mean pH  files

# Load BF output for 5 replicate Baypass runs
baypass.ph <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run", i, "_summary_yij_pij.out", sep=""))
    tmp[,rep:=i]
    return(tmp)
}

# ================================================================================== #

# Average  across replicate runs
# M_P is the mean of the posterior dist of the aij param, which is closely related to the freq of the reference allele
# M_Y is the posterior mean of allele counts)
baypass.ph.sum <- baypass.ph %>% group_by(POP, MRK) %>% 
    reframe(M_Y.mean = mean(M_Y), M_P.mean = mean(M_P))

# Join pops with Baypass output
baypass.ph.sum <- left_join(baypass.ph.sum, pops, by = "POP")

# ================================================================================== #

# Save output
save(baypass.ph.sum, "data/processed/baypass/abiotic/baypass.ph.sum.RData")
