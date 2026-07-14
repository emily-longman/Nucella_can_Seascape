# Nucella ABC analysis
## Part 1 - real data

# Clear memory
rm(list=ls())

# Stop exponential
options(scipen = 999)

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'lme4'))
library(data.table)
library(tidyverse)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(lme4)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load and format data

# Ecological variables
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Allele frequency data for top ph hits
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit and join with eco vars
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595")  %>%
  left_join(dplyr::select(ecovars, Site, Latitude, sim_eq)) %>%
  arrange(-Latitude)

# Make a pooled object with the raw data
pool.real <- new("pooldata",
                 npools=19, #### Rows = Number of pools
                 nsnp=1, ### Columns = Number of SNPs
                 refallele.readcount=as.matrix(t(topsnp$Count)),
                 readcoverage=as.matrix(t(topsnp$Cov)),
                 poolsizes=topsnp$nsnails * 2,
                 poolnames = topsnp$Site)

# ================================================================================== #

# Estimate the key statistics

# Hierach fst between N vs S
fst.phylogeo <- computeFST(pool.real,
                        method = "Anova",
                        struct = topsnp$shape, verbose = FALSE)
# "FST": estimate of genome-wide Fst over all the populations
# "FSG": estimate of genome-wide within-group differentiation (Fsg)
# "FGT": estimate of genome-wide between-group differentiation (Fgt)

# Raw correlation between mean pH and AF
rawcor = cor.test(~ ph_mean+AF, data = topsnp)

# Correlation between AF and AF - 1 for real data
rawcorAF = cor.test(topsnp$AF, topsnp$AF)

# ================================================================================== #

# Extract AFs of top SNP and format similar to sim output
AF_pool = data.frame(t(topsnp$AF))
names(AF_pool) = topsnp$sim_eq

# Format Data
real_data =
data.frame(
  Fsg = fst.phylogeo$snp.Fstats[1],
  Fgt = fst.phylogeo$snp.Fstats[2],
  Fst = fst.phylogeo$snp.Fstats[3],
  cor = rawcor$estimate,
  corAF = rawcorAF$estimate,
  fix1 = sum(topsnp$AF ==1),
  fix0 = sum(topsnp$AF ==0),
  poly = sum(topsnp$AF  > 0 & topsnp$AF  < 1),
  mean.AF = mean(topsnp$AF),
  AF_pool)

# ================================================================================== #

# Save output
save(real_data, file = "data/processed/SLiM/ph_ABC/real_data.Rdata")
