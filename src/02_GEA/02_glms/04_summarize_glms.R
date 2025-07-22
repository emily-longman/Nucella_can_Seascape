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
#install.packages(c('data.table', 'tidyverse', 'foreach'))
library(data.table)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Generate output directories

out_dir <- paste("data/processed/GEA/glms/glms_output_by_model")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load
load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")
#load("data/processed/GEA/glms/glms_chunk_analysis/GLM_100perm_Bio-Oracle_chunk_1.Rdata") # Load just one chunk to test

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Formatting

# Get names of enviro variables
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# Create SNP_id column
glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Test purposes
glm.model.collated <- glm.model.collated[which(glm.model.collated$chunk =="chunk_1.Rdata"),]

# Create separate dataframes for the real data the permutation data
real_data <- glm.model.collated %>% filter(perm == 0) #62,118
perm_data <- glm.model.collated %>% filter(perm > 0) #621,180

perm_data_sum <- foreach(i:length(enviro_vars_names), .combine="rbind")%do%{
    # Get variable name
    var <- enviro_vars_names[i]
    # Filter model output for just that model
    tmp <- perm_data %>% filter(variable == var)
    tmp_sum <- tmp %>% group_by(SNP_id) %>% 
            reframe(
            SNP_id=SNP_id,
            chr=chr,
            pos=pos,
            variable=variable,
            mean.AIC=mean(AIC),
            prop.perm.mu=mean(p_lrt),
            prop.perm.lci=quantile(p_lrt, .01),
            prop.perm.uci_0.95=quantile(p_lrt, .95),
            prop.perm.uci_0.99=quantile(p_lrt, .99),
            prop.perm.med=median(p_lrt)) %>% as.data.frame() %>% distinct()
}

# Join summary of perm data with real data
data_sum <- left_join(real_data, perm_data_sum, by = join_by(chr, pos, variable, SNP_id))


# ALTERNATE
glm_all_sum <- foreach(i:length(enviro_vars_names), .combine="rbind")%do%{
    # Get variable name
    var <- enviro_vars_names[i]
    # Filter model output for just that model
    tmp <- glm.model.collated %>% filter(variable == var)
    tmp_sum <- tmp %>% group_by(SNP_id, data) %>% 
            reframe(
            SNP_id=SNP_id,
            chr=chr,
            pos=pos,
            data=data,
            variable=variable,
            mean.AIC=mean(AIC),
            mu=mean(p_lrt),
            lci_0.01=quantile(p_lrt, .01),
            uci_0.95=quantile(p_lrt, .95),
            uci_0.99=quantile(p_lrt, .99),
            med=median(p_lrt)) %>% as.data.frame() %>% distinct()
}



# ================================================================================== #

# Create unique Chr number
chr.unique <- unique(glm_all_sum$chr)
glm_all_sum$chr.unique <- as.numeric(factor(glm_all_sum$chr, levels = chr.unique))

# For testing just extract one var
data_sum_all_temp <- glm_all_sum %>% filter(variable == "thetao_max")

# Graph - test
pdf("output/figures/GEA/glms/glm_perm_test.pdf", width = 8, height = 8)
ggplot(data_sum_all_temp, aes(x=chr.unique, y=-log(lci_0.01), color=data)) + geom_point() +
geom_hline(yintercept = -log10(0.01)) +
theme_bw()
dev.off()
