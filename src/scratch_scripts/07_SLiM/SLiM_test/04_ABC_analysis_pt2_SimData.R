# Nucella ABC analysis
## Part 2 - sim data

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'foreach', 'doParallel', 'lme4', 'abc'))
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

# Allele frequency data for top ph hits
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit and join with ecovars
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595") %>%
  left_join(dplyr::select(ecovars, Site, ph_mean))

# Load sim data
sim_Data <- get(load("data/processed/SLiM/ph_results.Rdata"))
names(sim_Data)[1:19] = paste("p", 0:18, sep ="")
# Reformat
sim_Data.melt <- sim_Data %>%
  reshape2::melt(id = c("repId","m","thresh","k_1","k_2",
                        "N","state","sim.cycle"),
                 variable.name = "sim_eq",
                 value.name = "AF_true")

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

# Create function that generates poolseq noise for sims, estimates key statistics, and formats data
# Note: must also load ecovars.sims, topsnp, sim_Data.melt

# Create function
process_sims = function(r, m, t, k1, k2, N){
  #r=1; m=0.001; t=7.96; k1=0; k2=0.001; N=2500

  # Extract data for specific parameter combos
  tmp <- sim_Data.melt %>%
    filter(repId==r,
           m==m,
           thresh==t,
           k_1==k1,
           k_2==k2,
           N==N)

  # Join with ecovars and top snp data
  tmp2 <- left_join(tmp, ecovars.sims, by = join_by(sim_eq)) %>%
          left_join(dplyr::select(topsnp, Site, nsnails, Cov), by = join_by(Site))

  ####

  # Step1-generate poolseq noise
  tmp.pool <- tmp2 %>%
    group_by(sim_eq) %>%
    mutate(SIM_AD = rbinom(1, Cov, AF_true)) %>%
    mutate(SIM_AF = SIM_AD/Cov)

  # Create pooled object
  pool.sim <- new("pooldata",
                   npools=19, #### Rows = Number of pools
                   nsnp=1, ### Columns = Number of SNPs
                   refallele.readcount=as.matrix(t(tmp.pool$SIM_AD)),
                   readcoverage=as.matrix(t(tmp.pool$Cov)),
                   poolsizes=tmp.pool$nsnails * 2,
                   poolnames = tmp.pool$Site )

  ####

  # Estimate the key statistics

  # Hierach fst
  fst.phylogeo <- computeFST (pool.sim,
                              method = "Anova",
                              struct = tmp.pool$`Demographic Cluster`, verbose = FALSE)

  #tmp.pool %<>%
  #  mutate(nEff =round((Cov*nsnails)/(Cov+nsnails-1)) ) %>%
  #  mutate(af_nEff:=round(SIM_AF*nEff)/nEff)
  #GLM = glmer(cbind(af_nEff*nEff, (1-af_nEff)*nEff) ~ ph_mean + (1 | `Demographic Cluster`),  
  #            data=tmp.pool,  family = binomial)
  #GLM_s = summary(GLM)

  # Raw correlation between mean pH and SIM AF
  rawcor = cor.test(~ ph_mean+SIM_AF, data =  tmp.pool)
  
  # Correlation between AF of the real data and AF of the sim data
  rawcorAF = cor.test(topsnp$AF, tmp.pool$SIM_AF)

  ## AFs themselves
  #AF_line = t(tmp.pool[,c("af_nEff")])
  #colnames(AF_line) = tmp.pool$Site
  
  # Extract AFs of top SNP and format similar to sim output
  AF_pool = data.frame(t(tmp.pool$AF_true))
  names(AF_pool) = tmp.pool$sim_eq

  ####

  # Format Data
  sim.data =
    data.frame(
      Fsg = fst.phylogeo$snp.Fstats[1],
      Fgt = fst.phylogeo$snp.Fstats[2],
      Fst = fst.phylogeo$snp.Fstats[3],
      #Beta = GLM_s$coefficients[2,1],
      cor = rawcor$estimate,
      corAF = rawcorAF$estimate,
      fix1 = sum(tmp.pool$SIM_AF ==1),
      fix0 = sum(tmp.pool$SIM_AF ==0),
      poly = sum(tmp.pool$SIM_AF  > 0 & topsnp$SIM_AF  < 1),
      mean.AF = mean(tmp.pool$SIM_AF),
      AF_pool
    ) 
  
  return(sim.data)
}

# ================================================================================== #

# Deploy function

# Extract parameters
repId = unique(sim_Data.melt$repId)
m = unique(sim_Data.melt$m)
thresh= unique(sim_Data.melt$thresh)
k_1 = unique(sim_Data.melt$k_1)
k_2 = unique(sim_Data.melt$k_2)
N = unique(sim_Data.melt$N)

#### This is a HIDEOUS LOOP!!
#### this must be optimized somehow ... but for now this will do...
results = 
foreach(r=repId,
        .errorhandling = "remove", .combine = "rbind")%do%{
  foreach(m=m,
          .errorhandling = "remove", .combine = "rbind")%do%{ 
    foreach(t=thresh,
            .errorhandling = "remove", .combine = "rbind")%do%{ 
      foreach(k1=k_1,
              .errorhandling = "remove", .combine = "rbind")%do%{ 
        foreach(k2=k_2,
                .errorhandling = "remove", .combine = "rbind")%do%{ 
          foreach(N=N,
                  .errorhandling = "remove", .combine = "rbind")%do%{ 
                    
   process_sims(r, m, t, k1, k2, N) %>%
                      mutate(
                        repId = r,
                        m = m,
                        thresh= t,
                        k_1 = k1,
                        k_2 = k2,
                        N = N
                      )
}}}}}}

# ================================================================================== #

# Save output
save(results, file = "data/processed/SLiM/ph_ABC/sim.results.Rdata")

#pdf("output/figures/SLiM/pH_sim_results_corAF.pdf", width = 5, height = 5)
#ggplot(results, aes(x=corAF)) + geom_density()
#dev.off()