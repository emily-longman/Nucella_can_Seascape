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

# Install and load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
w = as.numeric(args[1]) # Chunk: this is the chunk (1000 chunks each with 19 scaffolds in it)

# ================================================================================== #

# Load data

# Load txt file with scaffold names
scaffold.names.df <- read.csv(paste("data/processed/GEA/glms/scaffold.names", w, "txt", sep = "."), sep = " ", header=F)
scaffold.names <- scaffold.names.df$V1

# Load SNPs of interest (baypass POD outlier SNPs - 3,095 SNPs SNPs)
#baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs", header=T)
baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.01", header=T)
#baypass_POD_sig_SNPs <- read.table("data/processed/outlier_analyses/baypass/POD/baypass_POD_sig_SNPs_threshold_0.001", header=T)

# Load MARINe environmental data
# Note: data was re-arranged in excel to add my site names and do slight reformatting
marine_data <- read.csv("data/processed/GEA/enviro_data/MARINe/MARINe_data.csv", header=T)

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Load PCA data for demography
pca.df <- read.csv("data/processed/outlier_analyses/pca.csv")
colnames(pca.df)[1] <- "sampleId"

# ================================================================================== #

# Format data

# Rename "location" column as "sampleId"
names(marine_data)[names(marine_data) == "location"] <- "sampleId"

# Join bio-oracle dataframe and PCA dataframe
marine_data <- dplyr::left_join(marine_data, pca.df, by = "sampleId")

# Create SNP_id column for outlier SNP list
baypass_POD_sig_SNPs <- baypass_POD_sig_SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

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
anovaFun <- function(m1, m2) {
  ll1 <- as.numeric(logLik(m1))
  ll2 <- as.numeric(logLik(m2))
  parameter <- abs(attr(logLik(m1), "df") -  attr(logLik(m2), "df"))
  chisq <- -2*(ll1-ll2)
  1-pchisq(chisq, parameter)
    }

# ================================================================================== #
# ================================================================================== #

# Filter snp.dt for a given chunk - i.e., identify all of the sig SNPs in a given set of scaffold names
snp.dt %>% filter(snp.dt$chr %in% scaffold.names) -> data_chunk

# ================================================================================== #

glm.model.output =

  # For each SNP in a given chunk w extract allele freq then run model with bio-oracle data
  foreach(i=1:dim(data_chunk)[1], .combine = "rbind")%do%{
    
    # Reset filter
    seqResetFilter(genofile)

    ###############################################################

    # Calculate allele frequency for SNP i in chunk w 
    seqSetFilter(genofile, variant.id=data_chunk$variant.id[i], verbose = F)

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

    # Join with bio-oracle environmental data
    left_join(af_i_snp, marine_data, by ="sampleId") -> af_i_snp_enviro
    
    # Create long format data table with the enviro data in column labeled "value" and the specific variable identified in the column "column"
    af_i_snp_enviro %>% as_tibble %>% gather(key = "enviro_var", value = "value", `B.glandula_mean_perc_cov`:`P.ochraceus_mean_density_post_SSWD`) -> gathered_data
    
    # Get names of environmental variables
    unique(gathered_data$enviro_var) -> enviro_vars_names
    
    ###############################################################
  
    # Run model for each variable
    real_estimates =
      foreach(j=enviro_vars_names, .combine = "rbind", .errorhandling = "remove")%do%{
        
        # Extract data for 'j' environmental variable
        gathered_data %>% filter(enviro_var == j) -> inner.tmp
        
        # Remove rows when MARINe data is NA
        inner.tmp.sub <- inner.tmp[!is.na(inner.tmp$value),]

        # Model allele freq
        # Generate 3 models - a null model (t0), a model with just demography (t1.dem) and a model with demography and "j" enviro var (t1.dem.env)
        y <- inner.tmp.sub$af_nEff
        X.null <- model.matrix(~1, inner.tmp.sub)
        X.dem <- model.matrix(~PC1, inner.tmp.sub)
        X.dem.env <- model.matrix(~PC1+value, inner.tmp.sub)
        t0 <- fastglm(x=X.null, y=y, family=binomial(), weights=inner.tmp.sub$nEff, method=0)
        t1.dem <- fastglm(x=X.dem, y=y, family=binomial(), weights=inner.tmp.sub$nEff, method=0)
        t1.dem.env <- fastglm(x=X.dem.env, y=y, family=binomial(), weights=inner.tmp.sub$nEff, method=0)
        
        # Generate output table with t1.dem.env model information and model comparison info for each variable
        data.frame(
          chr = unique(inner.tmp.sub$chr),
          pos = unique(inner.tmp.sub$pos),
          variable = j,
          missing = seqMissing(genofile),
          perm = 0,
          data = "real",
          AIC_null = c(AIC(t0)),
          AIC_dem = c(AIC(t1.dem)),
          AIC_dem_env = c(AIC(t1.dem.env)),
          b_enviro = last(t1.dem.env$coef),
          se_enviro = last(t1.dem.env$se),
          p_lrt = anovaFun(t1.dem, t1.dem.env))
      } # End for enviro var

      ###############################################################
      
    # Permutations to generate null expectation of association between af and enviro var
    permutation_estimates =
      foreach(j=enviro_vars_names, .combine = "rbind", .errorhandling = "remove")%do%{
                
        # Extract data for 'j' environmental variable
        gathered_data %>% filter(enviro_var == j) -> inner.tmp.shuffle

          # Remove rows when MARINe data is NA
          inner.tmp.shuffle.sub <- inner.tmp.shuffle[!is.na(inner.tmp.shuffle$value),]
          
          # Do 100 permutations - 50 at a time to make things run faster (1:50 and 51:100)
          foreach(l=1:50, .combine = "rbind")%do%{
            set.seed(l)

            # Shuffle enviro data for 'j' enviro variable
            #inner.tmp.perm %>% mutate(shuffle_value = sample(value)) -> inner.tmp.shuffle
            inner.tmp.shuffle.sub$shuffle_value <- sample(inner.tmp.shuffle.sub$value)
            
            # Model allele freq
            # Generate 3 models - a null model (t0), a model with just demography (t1.dem) and a model with demography and "j" enviro var (t1.dem.env)
            y.perm <- inner.tmp.shuffle.sub$af_nEff
            X.null.perm <- model.matrix(~1, inner.tmp.shuffle.sub)
            X.dem.perm <- model.matrix(~PC1, inner.tmp.shuffle.sub)
            X.dem.env.perm <- model.matrix(~PC1+shuffle_value, inner.tmp.shuffle.sub)
            t0.perm <- fastglm(x=X.null.perm, y=y.perm, family=binomial(), weights=inner.tmp.shuffle.sub$nEff, method=0)
            t1.dem.perm <- fastglm(x=X.dem.perm, y=y.perm, family=binomial(), weights=inner.tmp.shuffle.sub$nEff, method=0)
            t1.dem.env.perm <- fastglm(x=X.dem.env.perm, y=y.perm, family=binomial(), weights=inner.tmp.shuffle.sub$nEff, method=0)
                  
                # Generate output table with t1.dem.env model information and model comparison info for each variable
                data.frame(
                  chr = unique(inner.tmp.shuffle.sub$chr),
                  pos = unique(inner.tmp.shuffle.sub$pos),
                  variable = j,
                  missing = seqMissing(genofile),
                  perm = l,
                  data = "permutation",
                  AIC_null = c(AIC(t0.perm)),
                  AIC_dem = c(AIC(t1.dem.perm)),
                  AIC_dem_env = c(AIC(t1.dem.env.perm)),
                  b_enviro = last(t1.dem.env.perm$coef),
                  se_enviro = last(t1.dem.env.perm$se),
                  p_lrt = anovaFun(t1.dem.perm, t1.dem.env.perm))

              } # End for perm
            } # End for enviro var
   
    # Combine real estimates and permutations
    all_data <- rbind(real_estimates, permutation_estimates)
    return(all_data)
  }

# ================================================================================== #

# Generate folders and save output

# Folder name for chunk c
folder_name <- paste("data/processed/GEA/glms/glms_chunk_analysis_marine")

# Save file for chunk w
file_name <- paste("GLM_50perm_marine_chunk_", w, sep = "")
save(glm.model.output, file = paste(folder_name, "/" , file_name, ".Rdata", sep = "") )

message("done")