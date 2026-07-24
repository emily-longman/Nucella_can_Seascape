# Merge SLiM output - ABC analysis sim data

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'foreach', 'doParallel', 'lme4', 'abc', 'stats', 'minpack.lm'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(foreach)
library(doParallel)
library(lme4)
library(abc)
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

# Allele frequency data for top ph hits
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit and join with ecovars
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595") %>%
  left_join(dplyr::select(ecovars, Site, ph_mean))

# ================================================================================== #

# Assigned ph values
env = c(
  8.023377561,
  8.016675575,
  8.013330447,
  8.006025013,
  7.987332623,
  7.926713508,
  7.957538759,
  7.940880291,
  7.937233951,
  7.933377415,
  7.923570354,
  7.941466321,
  7.941527015,
  7.946801475,
  7.963869108,
  7.96609182,
  7.948259634,
  7.946165852,
  7.950069857);

# Merge with ecovar
ecovars.sims <- cbind(ecovars, env) %>%
  mutate(sim_eq = paste("p", 0:18, sep =""))

# ================================================================================== #
# ================================================================================== #

# Function with more variables

# Create function
process_sims = function(repId, m, thresh, k, mag, N){
  #repId=1; m=0.005; thresh=7.985; k=0.04; mag=1; N=5000

  # Extract data for specific parameter combos
  tmp <- sim_Data_melt %>%
    filter(repId==repId,
           m==m,
           thresh==thresh,
           k==k,
           mag==mag,
           N==N)

  # Join with ecovars and top snp data
  tmp2 <- left_join(tmp, ecovars.sims, by = join_by(sim_eq)) %>%
          left_join(dplyr::select(topsnp, Site, nsnails, Cov), by = join_by(Site))

  #### Create poolobject and generate poolseq noise
  
  # Step1-generate poolseq noise
  tmp.pool <- tmp2 %>%
    group_by(sim_eq) %>%
    mutate(SIM_AD = rbinom(1, Cov, AF_true)) %>%
    mutate(SIM_AF = SIM_AD/Cov)
  
  # Create variable for AF groups
  tmp.pool$group <- c("L", "L", "L", "L", "M", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H")

  # Create pooled object
  pool.sim <- new("pooldata",
                   npools=19, #### Rows = Number of pools
                   nsnp=1, ### Columns = Number of SNPs
                   refallele.readcount=as.matrix(t(tmp.pool$SIM_AD)),
                   readcoverage=as.matrix(t(tmp.pool$Cov)),
                   poolsizes=tmp.pool$nsnails * 2,
                   poolnames = tmp.pool$Site )

  #### Estimate the key statistics

  # Hierach fst
  fst.phylogeo <- computeFST (pool.sim,
                              method = "Anova",
                              struct = tmp.pool$`Demographic Cluster`, verbose = FALSE)
  

  # Fst between pops with low and med AF
  pool.sim.L.M <- new("pooldata",
                   npools=5, #### Rows = Number of pools
                   nsnp=1, ### Columns = Number of SNPs
                   refallele.readcount=as.matrix(t(tmp.pool$SIM_AD[1:5])),
                   readcoverage=as.matrix(t(tmp.pool$Cov[1:5])),
                   poolsizes=tmp.pool$nsnails[1:5] * 2,
                   poolnames = tmp.pool$Site[1:5] )
  fst.L.M.group <- computeFST (pool.sim.H.M,
                              method = "Anova",
                              struct = tmp.pool$group[1:5], verbose = FALSE)
  # Fst between pops with med and high AF
  pool.sim.M.H <- new("pooldata",
                   npools=15, #### Rows = Number of pools
                   nsnp=1, ### Columns = Number of SNPs
                   refallele.readcount=as.matrix(t(tmp.pool$SIM_AD[5:19])),
                   readcoverage=as.matrix(t(tmp.pool$Cov[5:19])),
                   poolsizes=tmp.pool$nsnails[5:19] * 2,
                   poolnames = tmp.pool$Site[5:19])
  fst.M.H.group <- computeFST (pool.sim.M.H,
                              method = "Anova",
                              struct = tmp.pool$group[5:19], verbose = FALSE)
  # Fst between pops with low and high AF
  pool.sim.L.H <- new("pooldata",
                   npools=18, #### Rows = Number of pools
                   nsnp=1, ### Columns = Number of SNPs
                   refallele.readcount=as.matrix(t(tmp.pool$SIM_AD[-5])),
                   readcoverage=as.matrix(t(tmp.pool$Cov[-5])),
                   poolsizes=tmp.pool$nsnails[-5] * 2,
                   poolnames = tmp.pool$Site[-5])
  fst.L.H.group <- computeFST (pool.sim.L.H,
                              method = "Anova",
                              struct = tmp.pool$group[-5], verbose = FALSE)

  # Calculate the mean delta AF between groups
  mean.deltaAF.L.M <- mean(tmp.pool$SIM_AF[which(tmp.pool$group == "L")] - tmp.pool$SIM_AF[which(tmp.pool$group == "M")])
  mean.deltaAF.H.M <- mean(tmp.pool$SIM_AF[which(tmp.pool$group == "H")] - tmp.pool$SIM_AF[which(tmp.pool$group == "M")])

  # Raw correlation between mean pH and SIM AF
  rawcor.pearson = cor.test(~ ph_mean+SIM_AF, method = "pearson", data =  tmp.pool)
  rawcor.spearman = cor.test(~ ph_mean+SIM_AF, method = "spearman", data =  tmp.pool, exact = FALSE) #Goal is to cal correlation coef rho, not p-val

  # Correlation between AF of the real data and AF of the sim data
  rawcorAF.pearson = cor.test(topsnp$AF, tmp.pool$SIM_AF, method = "pearson")
  rawcorAF.spearman = cor.test(topsnp$AF, tmp.pool$SIM_AF, method = "spearman", exact = FALSE)

  # Fit sigmoid
  #tmp.pool_sub <- tmp.pool[,c("SIM_AF", "ph_mean")]
  #mod <- nlsLM(SIM_AF ~ SSlogis(ph_mean, Asym, xmid, scal), data = tmp.pool_sub)
  #mod_fit <- coef(mod)

  # Extract AFs of top SNP and format similar to sim output
  AF_pool = data.frame(t(tmp.pool$AF_true))
  names(AF_pool) = tmp.pool$sim_eq

  #### Format output

  # Format Data
  sim.data =
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
      #asym = mod_fit[1],
      #xmid = mod_fit[2],
      #scal = mod_fit[3],
      fix1 = sum(tmp.pool$SIM_AF ==1),
      fix0 = sum(tmp.pool$SIM_AF ==0),
      poly = sum(tmp.pool$SIM_AF > 0 & tmp.pool$SIM_AF < 1),
      mean.AF = mean(tmp.pool$SIM_AF),
      AF_pool)
  
  return(sim.data)
}

# ================================================================================== #

# Load data and run function

# Create list of file names for data
file_names = as.list(dir(path = 'data/processed/SLiM/ph_results_morevars5/', pattern = "phclineAFs.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/ph_results_morevars5/', x))))

# Read all the files and perform ABC
sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    #message(i)

    # Load file and add file name identifier
    tmp <- read.table(i) %>% 
            mutate(file_name = i) %>% 
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/ph_results_morevars5/phclineAFs.")) %>% 
            mutate(file_name = str_remove(file_name, pattern = ".txt"))
    # Rename pops
    names(tmp)[1:19] = paste("p", 0:18, sep ="")
    # Separate columns based on parameters
    sim_Data <- separate_wider_delim(tmp, cols = file_name, delim = "_", names = c("repId", "m", "thresh", "k", "mag", "N", "state", "sim.cycle"))

    # Reformat
    sim_Data_melt <- sim_Data %>%
            reshape2::melt(id = c("repId","m","thresh","k","mag","N","state","sim.cycle"), 
            variable.name = "sim_eq", value.name = "AF_true")
    
    # Extract parameters for file i
    repId_i = unique(sim_Data_melt$repId)
    m_i = unique(sim_Data_melt$m)
    thresh_i= unique(sim_Data_melt$thresh)
    k_i = unique(sim_Data_melt$k)
    mag_i = unique(sim_Data_melt$mag)
    N_i = unique(sim_Data_melt$N)

    # Run function
    process_sims(repId_i, m_i, thresh_i, k_i, mag_i, N_i) %>%
            mutate(repId = repId_i, m = m_i, thresh = thresh_i, k = k_i, mag = mag_i, N = N_i)

}

# ================================================================================== #

# Save output
#save(sim_data, file = "data/processed/SLiM/ph_ABC/sim_data_morevars2.Rdata")
#save(sim_data, file = "data/processed/SLiM/ph_ABC/sim_data_morevars3.Rdata")
save(sim_data, file = "data/processed/SLiM/ph_ABC/sim_data_morevars5.Rdata")