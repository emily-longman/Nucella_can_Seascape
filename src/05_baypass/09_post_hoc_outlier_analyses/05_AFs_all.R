# Calculate allele frequencies for entire SNP list

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

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# ================================================================================== #

# Calculate AF for every SNP

# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# ================================================================================== #

# Extract and manipulate count and coverage

# Extract read count data for SNPs
ref_count <- pooldata@refallele.readcount
ref_count %>% as.data.frame -> count
names(count) = c(pooldata@poolnames)
count$SNP_id <- snp.info$SNP_id
count.melt <- reshape2::melt(count, id = "SNP_id", variable.name = "Site", value.name = "Count")

# Extract coverage data for SNPs
coverage <- pooldata@readcoverage
coverage %>% as.data.frame -> cov
names(cov) = c(pooldata@poolnames)
cov$SNP_id <- snp.info$SNP_id
cov.melt <- reshape2::melt(cov, id = "SNP_id", variable.name = "Site", value.name = "Cov")

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage
# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata@poolnames)
# Add SNP_id
afs$SNP_id <- snp.info$SNP_id
# Change format
afs.melt <- reshape2::melt(afs, id = "SNP_id", variable.name = "Site", value.name = "AF")

# Join
afs.all <- left_join(count.melt, cov.melt)
afs.all <- left_join(afs.all, afs.melt)
afs.all <- left_join(snp.info, afs.all)

# Save
save(afs.all, file = "data/processed/outlier_analyses/afs.all.RData")
#load("data/processed/outlier_analyses/afs.all.RData")

# ================================================================================== #

# Summarize across sites to calculate global AFs
afs.all.ag <- afs.all %>% group_by(chr, pos, RefAllele, AltAllele, rs.id, SNP_id) %>%
    reframe(nSamps_poly=sum(AF>0 & AF<1, na.rm=T),
            nSamps_fixed=sum(AF==0 | AF==1, na.rm=T),
            global_af=sum(Count, na.rm=T)/sum(Cov, na.rm=T),
            mean_af=mean(AF, na.rm=T),
            poly_af=mean(AF[AF>0 & AF<1 & !is.na(AF)], na.rm=T))

# Save
save(afs.all.ag, file = "data/processed/outlier_analyses/afs.all.ag.RData")
#load("data/processed/outlier_analyses/afs.all.ag.RData")
