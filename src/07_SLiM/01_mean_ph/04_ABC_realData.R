# Nucella ABC analysis - real data

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'lme4', 'stats', 'minpack.lm'))
library(data.table)
library(tidyverse)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(lme4)
library(stats)
library(minpack.lm)

# ================================================================================== #

# Generate output directories

# Data directory
out_data_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

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
                        method = "Anova", struct = topsnp$shape, verbose = FALSE)
# "FST": estimate of genome-wide Fst over all the populations
# "FSG": estimate of genome-wide within-group differentiation (Fsg)
# "FGT": estimate of genome-wide between-group differentiation (Fgt)


fst.L.M.group <- computeFST (pool.sim.H.M,
                        method = "Anova", struct = tmp.pool$group[1:5], verbose = FALSE)

# Raw correlation between mean pH and AF
rawcor.pearson = cor.test(~ ph_mean+AF, method = "pearson", data = topsnp)
rawcor.spearman = cor.test(~ ph_mean+AF, method = "spearman", data = topsnp, exact = FALSE) #Goal is to cal correlation coef rho, not p-val

# Correlation between AF and AF - 1 for real data
rawcorAF.pearson = cor.test(topsnp$AF, topsnp$AF, method = "pearson",)
rawcorAF.spearman  = cor.test(topsnp$AF, topsnp$AF,  method = "spearman", exact = FALSE)

# Fit using self-starting parameters
topsnp_sub <- topsnp[,c(5,8)]
mod <- nlsLM(AF ~ SSlogis(ph_mean, Asym, xmid, scal), data = topsnp_sub)
mod_fit <- coef(mod)

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
  Fst.L.M = fst.L.M.group$Fst[1],
  Fst.M.H = fst.M.H.group$Fst[1],
  Fst.L.H = fst.L.H.group$Fst[1],
  deltaAF.L.M = mean.deltaAF.L.M, 
  deltaAF.H.M = mean.deltaAF.H.M,
  cor.pearson = rawcor.pearson$estimate,
  cor.spearman = rawcor.spearman$estimate,
  corAF.pearson = rawcorAF.pearson$estimate,
  corAF.spearman = rawcorAF.spearman$estimate,
  asym = mod_fit[1],
  xmid = mod_fit[2],
  scal = mod_fit[3],
  fix1 = sum(topsnp$AF ==1),
  fix0 = sum(topsnp$AF ==0),
  poly = sum(topsnp$AF  > 0 & topsnp$AF  < 1),
  mean.AF = mean(topsnp$AF),
  AF_pool)

# ================================================================================== #

# Save output
save(real_data, file = "data/processed/SLiM/ph_ABC/real_data.Rdata")
