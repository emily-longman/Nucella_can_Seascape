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
#load("data/processed/GEA/glms/glms_chunk_analysis_bio-oracle/glm.model.collated.subset1.Rdata") # 100 perm for SNPs in first half of the genome
# Load GLM biotic data
load("data/processed/GEA/glms/glms_chunk_analysis_marine/glm.model.collated.real.Rdata") # real data 
glm.biotic.real <- glm.model.collated
load("data/processed/GEA/glms/glms_chunk_analysis_marine/glm.model.collated.perm.1.25.Rdata") # 25 perm for SNPs
glm.biotic.perm.1.25 <- glm.model.collated
#load("data/processed/GEA/glms/glms_output/glm.model.collated.test.Rdata")
glm.biotic <- rbind(glm.biotic.real, glm.biotic.perm.1.25)

# Load bio-oracle environmental data
#bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# ================================================================================== #

# Format data

# Get names of enviro variables
#names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names
unique(glm.biotic.real$variable) -> enviro_vars_names

# Create SNP_id column
#glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
glm.biotic <- glm.biotic %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #


# Graph pval distribution - abiotic
#pdf("output/figures/GEA/glms/glm_pval_dist_abiotic.pdf", width = 8, height = 8)
#ggplot(glm.model.collated, aes(x=p_lrt, group=factor(perm), color=factor(perm))) + geom_density() +
#facet_wrap(~variable) + 
#scale_color_manual(values = c("red", rep("grey", 100))) +
#xlab("GLM P-values") + ylab("Number of SNPs") +
#theme_bw() + theme(legend.position = "none")
#dev.off()


# Graph pval distribution - biotic
#pdf("output/figures/GEA/glms/glm_pval_dist_biotic.pdf", width = 8, height = 8)
#ggplot(glm.model.collated, aes(x=p_lrt, group=factor(perm), color=factor(perm))) + geom_density() +
#facet_wrap(~variable) + 
#scale_color_manual(values = c("red", rep("grey", 25))) +
#xlab("GLM P-values") + ylab("Number of SNPs") +
#theme_bw() + theme(legend.position = "none")
#dev.off()

# ================================================================================== #

# Model enrichment

# Summarize the permutation data across the environmental variables
all_ratios <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
    # State variable name
    message(i)
    # Filter model output for just that model
    glm.model.collated %>% filter(variable == i) -> tmp_glm

    # Number of permutations
    n_perm = unique(tmp_glm$perm[which(tmp_glm$perm>0)])

    # Check if the environmental model beats the demographic model by comparing AIC
    tmp_glm <- tmp_glm %>% mutate(AICdiff = if_else(AIC_dem_env < AIC_dem, 1, 0)) 

    # Calculate relative rate (rr) of model enrichment (i.e., the number of SNPs where the env dem model was found as the better model)
    aic_sum <- tmp_glm %>% group_by(perm) %>% summarise(rr = sum(AICdiff == 1, na.rm = TRUE), .groups = "drop")

    # Separate by real and perm
    real_data <- aic_sum %>% filter(perm == 0) %>% rename(rr_real = rr)
    perm_data <- aic_sum %>% filter(perm > 0) %>% rename(rr_perm = rr)

    # rr of real data
    rr_real <- real_data$rr_real

    # Join the datasets and compare the rr between real and perm
    ratios <- perm_data %>%
    mutate(rr_ratio = log2(rr_real / rr_perm)) %>% mutate(variable = i)

}

# Write table
write.csv(all_ratios, "data/processed/GEA/glms/glms_summary/Biotic.rr.csv", row.names=FALSE)

# Summarize
summary <- all_ratios %>%
  group_by(variable) %>%
  drop_na() %>%
  summarise(mean = mean(rr_ratio), # Is the mean as simple as this or do I need to modify it?
            sd = sd(rr_ratio),
            n = n(),  # Number of observations
            se = sd / sqrt(n),  # Calculate standard error
            ci_low = mean - 1.96 * se,
            ci_high = mean + 1.96 * se,
            quantile_0.05 = quantile(rr_ratio, 0.05),
            quantile_0.95 = quantile(rr_ratio, 0.95),
            quantile_0.5 = quantile(rr_ratio, 0.5)) %>%
  ungroup()

# Write table
write.csv(summary, "data/processed/GEA/glms/glms_summary/Biotic.rr.sum.csv", row.names=FALSE)


# ================================================================================== #

# Graph summary
pdf("output/figures/GEA/glms/glm_biotic_rr_sum.pdf", width = 8, height = 8)
ggplot(summary, aes(x = reorder(variable, mean), y = mean)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd))+ 
  geom_point()+ 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  labs(title = "",
       x = "Group",
       y = "Value") +
  theme_minimal()+
  coord_flip()+
  theme_bw()
dev.off()


pdf("output/figures/GEA/glms/glm_biotic_rr_sum_CI.pdf", width = 8, height = 8)
ggplot(summary, aes(x = reorder(variable, mean), y = mean)) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high)) +  # Use 95% CI
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  labs(title = "Mean ± 95% CI Plot",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  coord_flip()
dev.off()

pdf("output/figures/GEA/glms/glm_biotic_rr_sum_quant_0.5.pdf", width = 8, height = 8)
ggplot(summary, aes(x = reorder(variable, quantile_0.5), y = quantile_0.5)) +
  geom_errorbar(aes(ymin = quantile_0.05, ymax = quantile_0.95)) +  # Use 95% CI
  geom_point() +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +  # Reference line at 1
  labs(title = "Mean ± 95% CI Plot",
       x = "Group",
       y = "Value") +
  theme_minimal() +
  coord_flip()
dev.off()

# ================================================================================== #
# ================================================================================== #

# Create separate dataframes for the real data and the permutation data
#real_data <- glm.model.collated %>% filter(perm == 0) #969705 (i.e., 9 * length(unique(real_data$SNP_id)))
#perm_data <- glm.model.collated %>% filter(perm > 0) #9697050 - makes sense bc did 10 permutations for these test runs

# Summarize the permutation data across the environmental variables
#perm_sum <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
#    # State variable name
#    message(i)
#    # Filter model output for just that model
#    perm_data %>% filter(variable == i) -> tmp_perm
#    tmp_perm %>% group_by(SNP_id) %>% 
#            reframe(
#            SNP_id=SNP_id,
#            chr=chr,
#            pos=pos,
#            data=data,
#            variable=variable,
#            perm_n=n(),
#            mean.AIC=mean(AIC),
#            perm.mu=mean(p_lrt),
#            perm.lci_0.01=quantile(p_lrt, .01),
#            perm.lci_0.05=quantile(p_lrt, .05),
#            perm.uci_0.95=quantile(p_lrt, .95),
#            perm.uci_0.99=quantile(p_lrt, .99),
#            perm.med=median(p_lrt)) %>% as.data.frame() %>% distinct()
#}

# Join real data and summary of permutations
#summary <- left_join(real_data, perm_sum, by = join_by(chr, pos, variable, SNP_id)) %>% select(!c(data.x, data.y))

# Does the real data beat permutations?
# T vs F
#summary <- summary %>% mutate(beat_perm = summary$p_lrt <= summary$perm.lci_0.01)
# 1 vs 0 
#summary <- summary %>% mutate(beat_perm_num = if_else(p_lrt <= perm.lci_0.01, 1, 0))



# ================================================================================== #

# Compare the environmental variables

# Count how many SNPs beat permutations for each env var
#env_sum <- foreach(i=enviro_vars_names, .combine="rbind")%do%{
#    message(i)

# Filter model output for just that model
#    summary %>% filter(variable == i) %>%
#            reframe(
#                variable=variable,
#                nSNPs=n(),
#                num_beat_perm = sum(beat_perm_num), 
#                binom.p = c(binom.test(num_beat_perm, nSNPs, p=0.05)$p.value)) %>% distinct()
#}
