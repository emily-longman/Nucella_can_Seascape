# Merge glms output

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
#install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)

# ================================================================================== #

# Generate output directories

out_dir <- paste("data/processed/GEA/glms/glms_output_by_model")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load
#load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")
load("data/processed/GEA/glms/glms_chunk_analysis/GLM_100perm_Bio-Oracle_chunk_1.Rdata") # Load just one chunk to test

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Get names of enviro variables
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# ================================================================================== #

# Create SNP_id column
glm.model.output <- glm.model.output %>% mutate(SNP_id = paste(chr, pos, sep = "_"))


o <- foreach(i:length(enviro_vars_names))%do%{
    # Get variable name
    var <- enviro_vars_names[i]

    # Filter model output for just that model
    tmp <- glm.model.output %>% filter(variable == var)

    # Remove af_nEff column
    tmp <- subset(tmp, select=-c(af_nEff))

    # Remove duplicated rows
    tmp <- tmp %>% distinct()

    o <- tmp %>% group_by(SNP_id, data) %>% summarize(
        SNP_id=SNP_id,
        variable=variable, 
        data=data,
        prop.perm.mu=median(p_lrt),
        prop.perm.lci=quantile(p_lrt, .01),
        prop.perm.uci=quantile(p_lrt, .99),
        prop.perm.med=median(p_lrt)
    )
}


o <- glm.model.output %>% group_by(SNP_id, variable) %>% mutate(
        prop=N/sum(N),
        prop.real=prop[perm==0], 
        totalN=totalN[perm==0],
        prop.perm.mu=median(prop[perm!=0]),
        prop.perm.lci=quantile(prop[perm!=0], .01),
        prop.perm.uci=quantile(prop[perm!=0], .99),
        prop.perm.med=median(prop[perm!=0]),
        prop.rr=median(log2(prop[perm==0]/prop[perm!=0])),
        prop.sd=sd(log2(prop[perm==0]/prop[perm!=0])))
            



### aggregate
  o2.ag <- o2[,list(prop.real=prop[perm==0], totalN=totalN[perm==0],
                    prop.perm.mu=median(prop[perm!=0]),
                    prop.perm.lci=quantile(prop[perm!=0], .01),
                    prop.perm.uci=quantile(prop[perm!=0], .99),
                    prop.perm.med=median(prop[perm!=0]),
                    prop.rr=median(log2(prop[perm==0]/prop[perm!=0])),
                    prop.sd=sd(log2(prop[perm==0]/prop[perm!=0]))),
                list(chr, inv=inv, mod, var, cluster, stage)]
  
  o2.ag[,rr:=prop.real/prop.perm.mu]
  o2.ag[,en:=(prop.real-prop.perm.mu)/prop.perm.mu]


  o2.ag[order(-prop.real)][]
  o2.ag[prop.real>(prop.perm.uci)][order(rr)]
  o2.ag[prop.rr-2*prop.sd>0]
  o2.ag[,p:=pnorm(0, prop.rr, prop.sd)]