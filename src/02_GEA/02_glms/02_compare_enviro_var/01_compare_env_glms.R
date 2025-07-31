# Summarize glms output and compare across the environmental variables

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories
out_dir <- paste("data/processed/GEA/glms/glms_summary")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load GLM data
load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")
#load("data/processed/GEA/glms/glms_output/glm.model.collated.test.Rdata")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Format data

# Get names of enviro variables
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# Create SNP_id column
glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Graph pval distribution
#pdf("output/figures/GEA/glms/glm_pval_dist.pdf", width = 8, height = 8)
#ggplot(glm.model.collated, aes(x=p_lrt, group=factor(perm), color=factor(perm))) + geom_density() +
#facet_wrap(~variable) + 
#scale_color_manual(values = c("red", rep("grey", 10))) +
#xlab("GLM P-values") + ylab("Number of SNPs") +
#theme_bw() + theme(legend.position = "none")
#dev.off()

# ================================================================================== #

# Create separate dataframes for the real data the permutation data
real_data <- glm.model.collated %>% filter(perm == 0) #969705 (i.e., 9 * length(unique(real_data$SNP_id)))
perm_data <- glm.model.collated %>% filter(perm > 0) #9697050 - makes sense bc did 10 permutations for these test runs

# Summarize the permutation data across the environmental variables
perm_sum <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
    # State variable name
    message(i)
    # Filter model output for just that model
    perm_data %>% filter(variable == i) -> tmp_perm
    tmp_perm %>% group_by(SNP_id) %>% 
            reframe(
            SNP_id=SNP_id,
            chr=chr,
            pos=pos,
            data=data,
            variable=variable,
            perm_n=n(),
            mean.AIC=mean(AIC),
            perm.mu=mean(p_lrt),
            perm.lci_0.01=quantile(p_lrt, .01),
            perm.lci_0.05=quantile(p_lrt, .05),
            perm.uci_0.95=quantile(p_lrt, .95),
            perm.uci_0.99=quantile(p_lrt, .99),
            perm.med=median(p_lrt)) %>% as.data.frame() %>% distinct()
}

# Join real data and summary of permutations
summary <- left_join(real_data, perm_sum, by = join_by(chr, pos, variable, SNP_id)) %>% select(!c(data.x, data.y))

# Does the real data beat permutations?
# T vs F
summary <- summary %>% mutate(beat_perm = summary$p_lrt <= summary$perm.lci_0.01)
# 1 vs 0 
summary <- summary %>% mutate(beat_perm_num = if_else(p_lrt <= perm.lci_0.01, 1, 0))

# Check to see how many beat perm
length(which(summary$beat_perm))

save(summary, file="data/processed/GEA/glms/glms_summary/summary.RData")
load("data/processed/GEA/glms/glms_summary/summary.RData")

# The dimensions of this table seem off - it should be the same size as real_data (i.e., 969705), but instead its 134077212

# ================================================================================== #

# Compare the environmental variables

# Count how many SNPs beat permutations for each env var
env_sum <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
    message(i)

# Filter model output for just that model
    summary %>% filter(variable == i) %>%
            reframe(
                variable=variable,
                nSNPs=n(),
                num_beat_perm = sum(beat_perm_num), 
                binom.p = c(binom.test(num_beat_perm, nSNPs, p=0.05)$p.value)) %>% distinct()
}

# ================================================================================== #

# Save summary
save(env_sum, file="data/processed/GEA/glms/glms_summary/env_sum.RData")
load("data/processed/GEA/glms/glms_summary/env_sum.RData")

# Problem - all are producing the same binomial p value

# ================================================================================== #

# Graph

# Create unique Chr number
#chr.unique <- unique(glm_all_sum$chr)
#glm_all_sum$chr.unique <- as.numeric(factor(glm_all_sum$chr, levels = chr.unique))

# For testing just extract one var
#data_sum_all_temp <- glm_all_sum %>% filter(variable == "thetao_max")

# Graph - test
#pdf("output/figures/GEA/glms/glm_perm_test.pdf", width = 8, height = 8)
#ggplot(data_sum_all_temp, aes(x=chr.unique, y=-log(mu), color=data)) + geom_point() +
#geom_hline(yintercept = -log10(0.01)) +
#theme_bw()
#dev.off()

# ================================================================================== #



