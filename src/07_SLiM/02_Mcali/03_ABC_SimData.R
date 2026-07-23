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
out_data_dir <- paste("data/processed/SLiM/Mcali_ABC")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# ================================================================================== #

# Load and format data

# Ecological variables
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Allele frequency data for top ph hits
afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T) %>% mutate(nsnails = 20)

# Order by site/latitude
afs.Mcali %<>% arrange(desc(latitude))

# ================================================================================== #

# Assigned Mcali values
env = c(
1.801197107,
1.758352318,
1.715225361,
1.677236687,
1.734608112,
1.844818925,
1.735326540,
1.615483995,
1.961925234,
1.738814502,
1.517563421,
1.487237983,
1.567637305,
1.472519437,
1.595566753,
2.206475463,
1.520675474,
1.431920910);

# Merge with ecovar
afs.Mcali.sims.tmp <- cbind(afs.Mcali, env)

# Join
afs.Mcali.sims <- left_join(afs.Mcali.sims.tmp, dplyr::select(ecovars, "Site", "sim_eq", "Demographic Cluster"))

# ================================================================================== #
# ================================================================================== #


# Function with more variables

# Create function
process_sims = function(repId, m, thresh, k, mag, N){
  #repId=1; m=0.005; thresh=1.7; k=0.2; mag=1; N=5000

  # Extract data for specific parameter combos
  tmp <- sim_Data_melt %>%
    filter(repId==repId,
           m==m,
           thresh==thresh,
           k==k,
           mag==mag,
           N==N)

  # Join with ecovars and top snp data
  tmp2 <- left_join(tmp, afs.Mcali.sims, by = join_by(sim_eq))

  #### Create poolobject and generate poolseq noise
  
  # Step1-generate poolseq noise
  tmp.pool <- tmp2 %>%
    group_by(sim_eq) %>%
    mutate(SIM_AD = rbinom(1, Cov, AF_true)) %>%
    mutate(SIM_AF = SIM_AD/Cov)

  # Create pooled object
  pool.sim <- new("pooldata",
                   npools=18, #### Rows = Number of pools
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

  # Raw correlation between shell thk and SIM AF
  rawcor.pearson = cor.test(~ mean_integrated_thk+SIM_AF, method = "pearson", data =  tmp.pool)
  rawcor.spearman = cor.test(~ mean_integrated_thk+SIM_AF, method = "spearman", data =  tmp.pool, exact = FALSE) #Goal is to cal correlation coef rho, not p-val
  
  # Correlation between AF of the real data and AF of the sim data
  rawcorAF.pearson = cor.test(afs.Mcali.sims$AF, tmp.pool$SIM_AF, method = "pearson")
  rawcorAF.spearman = cor.test(afs.Mcali.sims$AF, tmp.pool$SIM_AF, method = "spearman", exact = FALSE)

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
file_names = as.list(dir(path = 'data/processed/SLiM/Mcali_results_v3/', pattern = "McaliclineAFs.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/SLiM/Mcali_results_v3/', x))))

# Read all the files and perform ABC
sim_data <- foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    #message(i)

    # Load file and add file name identifier
    tmp <- read.table(i) %>% 
            mutate(file_name = i) %>% 
            mutate(file_name = str_remove(file_name, pattern = "data/processed/SLiM/Mcali_results_v3/McaliclineAFs.")) %>% 
            mutate(file_name = str_remove(file_name, pattern = ".txt"))
    # Rename pops
    names(tmp)[1:18] = paste("p", c(0,1,2,4:18), sep ="")
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
save(sim_data, file = "data/processed/SLiM/Mcali_ABC/sim_data_v3.Rdata")