# Graph glm results

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

# Get names of enviro variables
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)

# Extract just names
names(bio_oracle_sites_2010)[4:12] -> enviro_vars_names

# ================================================================================== #

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
