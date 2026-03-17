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
out_dir <- paste("data/processed/GEA/glms/glms_summary/model_enrichment")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load csv files

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_summary/', pattern = "Vars_rr_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_summary/', x))))

# Read all the files and add a column with the chunk
all_ratios <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    ratios = read.csv(w, header=T)
}

# Save all ratios
write.csv(all_ratios, "data/processed/GEA/glms/glms_summary/All_vars_rr.csv", row.names=FALSE)

# ================================================================================== #

# Get names of variables
abiotic_vars <- var_full_names$variable[which(var_full_names$group=="Abiotic")]
biotic_vars <- var_full_names$variable[which(var_full_names$group=="Biotic")]

# ================================================================================== #

# Summarize for abiotic vars
summary_abiotic <- all_ratios %>%
  filter(variable %in% abiotic_vars) %>%
  group_by(variable) %>%
  drop_na() %>%
  summarise(n = n(),  # Number of observations,
            mean_rr_log2 = mean(rr_ratio_log2),
            sd_rr_log2 = sd(rr_ratio_log2),
            se_rr_log2 = sd_rr_log2 / sqrt(n),  # Calculate standard error
            ci_low = mean_rr_log2 - 1.96 * se_rr_log2,
            ci_high = mean_rr_log2 + 1.96 * se_rr_log2,
            quantile_0.05 = quantile(rr_ratio_log2, 0.05),
            quantile_0.95 = quantile(rr_ratio_log2, 0.95),
            quantile_0.5 = quantile(rr_ratio_log2, 0.5)) %>%
  ungroup()
# Join with full names
summary_abiotic <- left_join(summary_abiotic, var_full_names, by="variable")
# Write table
write.csv(summary_abiotic, "data/processed/GEA/glms/glms_summary/Abiotic_vars_rr_sum.csv", row.names=FALSE)

# Summarize for abiotic vars
summary_biotic <- all_ratios %>%
  filter(variable %in% biotic_vars) %>%
  group_by(variable) %>%
  drop_na() %>%
  summarise(n = n(),  # Number of observations,
            mean_rr_log2 = mean(rr_ratio_log2),
            sd_rr_log2 = sd(rr_ratio_log2),
            se_rr_log2 = sd_rr_log2 / sqrt(n),  # Calculate standard error
            ci_low = mean_rr_log2 - 1.96 * se_rr_log2,
            ci_high = mean_rr_log2 + 1.96 * se_rr_log2,
            quantile_0.05 = quantile(rr_ratio_log2, 0.05),
            quantile_0.95 = quantile(rr_ratio_log2, 0.95),
            quantile_0.5 = quantile(rr_ratio_log2, 0.5)) %>%
  ungroup()
# Join with full names
summary_biotic <- left_join(summary_biotic, var_full_names, by="variable")
# Write table
write.csv(summary_biotic, "data/processed/GEA/glms/glms_summary/Biotic_vars_rr_sum.csv", row.names=FALSE)

# ================================================================================== #

# Graph summary

# Graph relative rate of model enrichment for abiotic - mean and 2*sd
pdf("output/figures/GEA/glms/model_enrichment/GLM_Abiotic_rr_sum_2SD.pdf", width = 6.75, height = 8)
ggplot(summary_abiotic, aes(x = reorder(variable_full_name, mean_rr_log2), y = mean_rr_log2)) + 
  geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2))+ 
  geom_point()+ 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  ylim(-0.65, 0.65) +
  labs(x = "",
       y = "Log2(Relative rate of model enrichment)") +
  theme_minimal()+
  coord_flip()+
  theme_bw(base_size=18)
dev.off()

# Graph relative rate of model enrichment for biotic - mean and 2*sd
pdf("output/figures/GEA/glms/model_enrichment/GLM_Biotic_rr_sum_2SD.pdf", width = 8, height = 8)
ggplot(summary_biotic, aes(x = reorder(variable_full_name, mean_rr_log2), y = mean_rr_log2)) + 
  geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2))+ 
  geom_point()+ 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
  ylim(-0.65, 0.65) +
  labs(x = "",
       y = "Log2(Relative rate of model enrichment)") +
  theme_minimal()+
  coord_flip()+
  theme_bw(base_size=18)
dev.off()

# ================================================================================== #

# Create group based on if beat 2sd or 1sd
summary_abiotic <- summary_abiotic %>% 
  mutate(col = case_when(mean_rr_log2-2*sd_rr_log2 > 0 ~ "blue",
                      mean_rr_log2-2*sd_rr_log2 < 0 & mean_rr_log2-sd_rr_log2 > 0 ~ "steelblue3",
                      mean_rr_log2-sd_rr_log2 < 0 ~ "gray56"))


# Graph relative rate of model enrichment for biotic - mean and 2*sd
pdf("output/figures/GEA/glms/model_enrichment/GLM_Abiotic_rr_color.pdf", width = 14, height = 5)
ggplot(summary_abiotic, aes(x = mean_rr_log2, y = reorder(variable_full_name, mean_rr_log2), color = col)) + 
  #geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2))+ 
  geom_point(size=7)+ 
  scale_color_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  xlim(-0.365, 0.365) + 
  labs(x = "Log2(Relative rate \nof model enrichment)",
       y = "") +
  theme_minimal()+
  coord_flip()+
  theme_bw(base_size=30) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
dev.off()
pdf("output/figures/GEA/glms/model_enrichment/GLM_Abiotic_rr_color_altsize.pdf", width = 9, height = 5)
ggplot(summary_abiotic, aes(x = mean_rr_log2, y = reorder(variable_full_name, mean_rr_log2), color = col)) + 
  #geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2))+ 
  geom_point(size=7)+ 
  scale_color_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  xlim(-0.365, 0.365) + 
  labs(x = "Log2(Relative rate \nof model enrichment)",
       y = "") +
  theme_minimal()+
  coord_flip()+
  theme_bw(base_size=30) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
dev.off()

# Create group based on if beat 2sd or 1sd
summary_biotic <- summary_biotic %>% 
  mutate(col = case_when(mean_rr_log2-2*sd_rr_log2 > 0 ~ "blue",
                      mean_rr_log2-2*sd_rr_log2 < 0 & mean_rr_log2-sd_rr_log2 > 0 ~ "steelblue3",
                      mean_rr_log2-sd_rr_log2 < 0 ~ "gray56"))

# Graph relative rate of model enrichment for biotic - mean and 2*sd
pdf("output/figures/GEA/glms/model_enrichment/GLM_Biotic_rr_color.pdf", width = 14, height = 5)
ggplot(summary_biotic, aes(x = mean_rr_log2, y = reorder(variable_full_name, mean_rr_log2), color = col)) + 
  #geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2))+ 
  geom_point(size=7)+ 
  scale_color_identity() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +  
  xlim(-0.365, 0.365) + 
  labs(x = "Log2(Relative rate \nof model enrichment)",
       y = "") +
  theme_minimal()+
  coord_flip()+
  theme_bw(base_size=30) + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())
dev.off()

# Graph relative rate of model enrichment - mean and sd
#pdf("output/figures/GEA/glms/GLM_rr_sum_SD.pdf", width = 8, height = 8)
#ggplot(summary, aes(x = reorder(variable, mean_rr_log2), y = mean_rr_log2)) + 
#  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2))+ 
#  geom_point()+ 
#  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  
#  labs(title = "Mean ± SD",
#       x = "",
#       y = "Log2(Relative rate of model enrichment)") +
#  theme_minimal()+
#  coord_flip()+
#  theme_bw(base_size=20)
#dev.off()

# Graph relative rate of model enrichment  - mean and 95% CI
#pdf("output/figures/GEA/glms/GLM_rr_sum_CI.pdf", width = 8, height = 8)
#ggplot(summary, aes(x = reorder(variable, mean_rr_log2), y = mean_rr_log2)) +
#  geom_errorbar(aes(ymin = ci_low, ymax = ci_high)) +  # Use 95% CI
#  geom_point() +
#  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
#  labs(title = "Mean ± 95% CI Plot",
#       x = "",
#       y = "Log2(Relative rate of model enrichment)") +
#  theme_minimal() +
#  coord_flip() + 
#  theme_bw(base_size=20)
#dev.off()