# Generalized linear models to assess relationship between SNPs and environmental variables

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
# install.packages(c('data.table', 'tidyverse', 'foreach', 'poolfstat', 'magrittr', 'reshape2', 'broom', 'stats', 'fastglm'))
library(data.table)
library(tidyverse)
library(foreach)
library(poolfstat)
library(magrittr)
library(reshape2)
library(broom)
library(stats)
library(fastglm)
library(poolfstat)

# Install and load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
w = as.numeric(args[1]) # Chunk: this is the chunk (1000 chunks each with 19 scaffolds in it)

# Prevent scientific notation
options(scipen=999)

# ================================================================================== #

# Load data

# Load txt file with scaffold names
scaffold.names.df <- read.csv(paste("data/processed/GEA/glms/scaffold.names", w, "txt", sep = "."), sep = " ", header=F)
scaffold.names <- scaffold.names.df$V1

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Load M californianus shell thickness data
Mcalifornianus_data <- read.csv("data/processed/GEA/enviro_data/Mcali_thk/Mcalifornianus_data_clean_18pop.csv", header=T)

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Load PCA data for demography
pca.df <- read.csv("data/processed/outlier_analyses/pca.csv")
colnames(pca.df)[1] <- "sampleId"
  
# ================================================================================== #

# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Extract SNP info for all SNPs
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  pooldata.snp.info

# Rename columns
names(pooldata.snp.info)[1:2] = c("chr","pos")

# Make snp_id column
pooldata.snp.info <- pooldata.snp.info %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Format data

# Rename "location" column as "sampleId"
names(Mcalifornianus_data)[names(Mcalifornianus_data) == "Site.Code"] <- "sampleId"
names(bio_oracle_sites_2010)[names(bio_oracle_sites_2010) == "location"] <- "sampleId"

# Extract just focal ecological vars
Mcalithk <- Mcalifornianus_data[,c(1,7)]
ph <- bio_oracle_sites_2010[,c(11,13)]

# Join datasets with PCA dataframe
ecological_data <- dplyr::left_join(Mcalithk, pca.df, by = "sampleId")
ecological_data <- dplyr::left_join(ph, ecological_data, by = "sampleId")

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Define anova function for comparing glms
#anovaFun <- function(m1, m2) {
#  ll1 <- as.numeric(logLik(m1))
#  ll2 <- as.numeric(logLik(m2))
#  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
#  chisq <- -2*(ll1-ll2)
#  1-pchisq(chisq, parameter)
#    }

# ================================================================================== #
# ================================================================================== #

# Filter snp.dt for a given chunk - i.e., identify all of the sig SNPs in a given set of scaffold names
snp.dt %>% filter(snp.dt$chr %in% scaffold.names) -> data_chunk

# Filter GLM to only sites in pooldata snp.info
data_chunk_filt <- data_chunk %>% filter(SNP_id %in% pooldata.snp.info$SNP_id)

# ================================================================================== #

glm.model.output =

  # For each SNP in a given chunk w extract allele freq then run model with bio-oracle data
  foreach(i=1:dim(data_chunk_filt)[1], .combine = "rbind")%do%{
    
    # Reset filter
    seqResetFilter(genofile, verbose = F)

    ###############################################################

    # Calculate allele frequency for SNP i in chunk w 
    seqSetFilter(genofile, variant.id=data_chunk_filt$variant.id[i], verbose = F)

    # Extract allele depth ('ad') of alternate allele for SNP i
    ad_i <- seqGetData(genofile, "annotation/format/AD") %>% .$data %>% .[,2]
    # Extract total depth ('dp') for SNP i
    dp_i <- seqGetData(genofile, "annotation/format/DP")[,1]

    # Create allele freq (af) data table with ad, dp, af, sample ID and variant id for SNP i
    af_i <- data.table(ad=ad_i, dp=dp_i, af=ad_i/dp_i,
    sampleId=seqGetData(genofile, "sample.id"),
    variant.id=rep(seqGetData(genofile, "variant.id"), each=length(ad_i)))

    # Merge allele freq table af_i and snp.dt
    af_i_snp <- merge(af_i, snp.dt, by="variant.id")

    # Sample size of each pool
    nSnail=20
    # Calculate the mean effective coverage ('nEff') (note: each pool consists of 20 dogwhelks)
    af_i_snp[,nEff:=round((dp*2*nSnail)/(2*nSnail+dp-1))]
    # Calculate the effective allele freq
    af_i_snp[,af_nEff:=round(af*nEff)/nEff]

    ###############################################################

    # Remove the two populations that are not in the MARINe database from af_i_snp
    af_i_snp <- af_i_snp[-which(af_i_snp$sampleId == "ARA" ),]

    ###############################################################

    # Join with environmental data
    left_join(af_i_snp, ecological_data, by ="sampleId") -> af_i_snp_enviro
      
    ###############################################################
  
    # Run model
        
        # Model allele freq
        # Generate 3 models - a null model (t0), a model with just demography (t1.dem) and a model with demography and "j" enviro var (t1.dem.env)
        y <- af_i_snp_enviro$af_nEff
        X.null <- model.matrix(~1, af_i_snp_enviro)
        X.dem <- model.matrix(~PC1, af_i_snp_enviro)
        X.dem.abiotic <- model.matrix(~PC1+ph_mean, af_i_snp_enviro)
        X.dem.biotic <- model.matrix(~PC1+mean_integrated_thk, af_i_snp_enviro)
        X.dem.both <- model.matrix(~PC1+ph_mean+mean_integrated_thk, af_i_snp_enviro)
        X.dem.both.int <- model.matrix(~PC1+ph_mean*mean_integrated_thk, af_i_snp_enviro)
        t0 <- fastglm(x=X.null, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        t1.dem <- fastglm(x=X.dem, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        t1.dem.abiotic <- fastglm(x=X.dem.abiotic, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        t1.dem.biotic <- fastglm(x=X.dem.biotic, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        t1.dem.both <- fastglm(x=X.dem.both, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        t1.dem.both.int <- fastglm(x=X.dem.both.int, y=y, family=binomial(), weights=af_i_snp_enviro$nEff, method=0)
        
        # Generate output table with t1.dem.env model information and model comparison info for each variable
        data.frame(
          chr = unique(af_i_snp_enviro$chr),
          pos = unique(af_i_snp_enviro$pos),
          SNP_id = unique(af_i_snp_enviro$SNP_id),
          variable = "abiotic_biotic",
          missing = seqMissing(genofile),
          perm = 0,
          data = "real",
          AIC_null = c(AIC(t0)),
          AIC_dem = c(AIC(t1.dem)),
          AIC_dem_abiotic = c(AIC(t1.dem.abiotic)),
          AIC_dem_biotic = c(AIC(t1.dem.biotic)),
          AIC_dem_both = c(AIC(t1.dem.both)),
          AIC_dem_both_int = c(AIC(t1.dem.both.int)),
          b_int = last(t1.dem.both.int$coef),
          se_int = last(t1.dem.both.int$se),
          p_int = summary(t1.dem.both.int)$coefficients[5,4])

      ###############################################################
 
  }

# ================================================================================== #

# Generate folders and save output

# Folder name
folder_name <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/real")
if (!dir.exists(folder_name)) {dir.create(folder_name)}

# Save file for chunk w
file_name <- paste("GLM_chunk_", w, sep = "")
save(glm.model.output, file = paste(folder_name, "/" , file_name, ".Rdata", sep = "") )

message("done")