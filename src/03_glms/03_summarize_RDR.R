# Summarize relative deviance reduction

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
#install.packages(c('data.table', 'tidyverse', 'dplyr', 'foreach', 'magrittr', 'RColorBrewer', 'ggtext))
library(data.table)
library(tidyverse)
library(dplyr)
library(foreach)
library(magrittr)
library(RColorBrewer)
library(ggtext)

# ================================================================================== #

# Generate output directories
out_fig_dir <- paste("output/figures/glms")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load seascape var full names
var_full_names <- read.csv("guide_files/Seascape_vars_full_names.csv", header=T)

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = "data/processed/glms/RDR", pattern = "glm_RDR_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/glms/RDR/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
RDR =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))

    # Return
    o
}

# Check structure
str(RDR)

# Change f2 to RDR
names(RDR)[names(RDR) == "f2"] <- "RDR"

# ================================================================================== #

# List of variables
vars <- unique(RDR$variable)

# Calc summary stats
effectsize_sum <- foreach(i = vars, .combine = "rbind")%do%{
    # Extract for just one var
    tmp <- RDR %>% filter(variable == i)
    # Summarize
    q <- quantile(tmp$RDR, probs = c(0.9, 0.95, 0.99, 0.995, 0.999), na.rm = T)
    q_log_pos <- quantile(log(tmp$RDR[tmp$RDR > 0]), probs = c(0.9, 0.95, 0.99, 0.995, 0.999), na.rm = T)
    
    # Make table
    data.frame(
        variable = i,
        n_tot = dim(tmp)[1],
        sum_RDR_neg = sum(tmp$RDR < 0),
        sum_RDR_0 = sum(tmp$RDR == 0),
        n_sig_0.05 = sum(tmp$p_lrt < 0.05),
        n_sig_0.01 = sum(tmp$p_lrt < 0.01),
        n_sig_0.001 = sum(tmp$p_lrt < 0.01/dim(tmp)[1]),
        mean_RDR_0.05 = mean(tmp$RDR[which(tmp$p_lrt < 0.05)]),
        median_RDR_0.05 = median(tmp$RDR[which(tmp$p_lrt < 0.05)]),
        mean_RDR_0.01 = mean(tmp$RDR[which(tmp$p_lrt < 0.01)]),
        median_RDR_0.01 = median(tmp$RDR[which(tmp$p_lrt < 0.01)]),
        mean_RDR_0.001 = mean(tmp$RDR[which(tmp$p_lrt < 0.001)]),
        median_RDR_0.001 = median(tmp$RDR[which(tmp$p_lrt < 0.001)]),
        mean_RDR =  mean(tmp$RDR),
        median_RDR =  median(tmp$RDR),
        mean_pos_RDR =  mean(tmp$RDR[tmp$RDR > 0]),
        median_pos_RDR =  median(tmp$RDR[tmp$RDR > 0]),
        RDR_quant_90 = q[1],
        RDR_quant_95 = q[2],
        RDR_quant_99 = q[3],
        RDR_quant_99.5 = q[4],
        RDR_quant_99.9 = q[5],
        log_pos_RDR_quant_90 = q_log_pos[1],
        log_pos_RDR_quant_95 = q_log_pos[2],
        log_pos_RDR_quant_99 = q_log_pos[3],
        log_pos_RDR_quant_99.5 = q_log_pos[4],
        log_pos_RDR_quant_99.9 = q_log_pos[5]
    )
}

# Join with full names
effectsize_sum <- left_join(effectsize_sum, var_full_names, by="variable")

# ================================================================================== #

# Graph effect size results

# Graph RDR scores for 95%
pdf("output/figures/glms/glm_RDR_95perc_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = RDR_quant_95, y = reorder(variable_full_name_edit, RDR_quant_95))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("95% Percentile " * RDR), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()
# Graph RDR scores for 99.9%
pdf("output/figures/glms/glm_RDR_99.9perc_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = RDR_quant_99.9, y = reorder(variable_full_name_edit, RDR_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile " * RDR), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph log of pos RDR scores for 95%
pdf("output/figures/glms/glm_log_pos_RDR_95perc_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = log_pos_RDR_quant_95, y = reorder(variable_full_name_edit, log_pos_RDR_quant_95))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("95% Percentile (log" * RDR *")"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()# Graph log pos RDR scores for 99.9%
pdf("output/figures/glms/glm_log_pos_RDR_99.9perc_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = log_pos_RDR_quant_99.9, y = reorder(variable_full_name_edit, log_pos_RDR_quant_99.9))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("99.9% Percentile log( " * RDR *")"), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# Graph median RDR
pdf("output/figures/glms/glm_median_pos_RDR_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = median_pos_RDR, y = reorder(variable_full_name_edit, median_pos_RDR))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Median pos" * RDR), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()
# Graph mean RDR
pdf("output/figures/glms/glm_mean_pos_RDR_all.pdf", width = 9, height = 12)
ggplot(effectsize_sum, aes(x = mean_pos_RDR, y = reorder(variable_full_name_edit, mean_pos_RDR))) + 
  geom_point(size = 8) + 
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  labs(x = expression("Mean pos" * RDR), y = "") +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_text(size = 22))
dev.off()

# ================================================================================== #

# Join with relative rate of enrichment

# Add model enrichment
summary <- read.csv("data/processed/glms/glms_summary/All_vars_rr_sum.csv", header = T)
summary_abiotic <- read.csv("data/processed/glms/glms_summary/Abiotic_vars_rr_sum.csv", header = T)
summary_biotic <- read.csv("data/processed/glms/glms_summary/Biotic_vars_rr_sum.csv", header = T)

# Join
summary_F <- left_join(summary, effectsize_sum)

# For graphing make species names italics
summary_F$variable_full_name_italics <- c("*B. glandula* (hmean)", "*B. glandula* (mean)", "*C. dalli* (hmean)", "*C. dalli* (mean)",
"*M. californianus* (hmean)", "*M. californianus* (mean)", "*M. trossulus* (hmean)", "*M. trossulus* (mean)", 
"*N. canaliculata* (hmean)", "*N. canaliculata* (mean)","*N. ostrina* (hmean)", "*N. ostrina* (mean)", 
"*P. ochraceus* (hmean)", "*P. ochraceus* (mean)", "Chlorophyll (mean)", "STI (mean)", "Cross-sectional thk (mean)", "Max thk (mean)", "Min thk (mean)",
"O2 (mean)", "pH (mean)", "pH (min)", "Salinity (mean)", "SST (max)", "SST (mean)", "SST (min)", "SST (range)")

# Graph relative rate of model enrichment for all vars - mean and 1*sd
pdf("output/figures/glms/GLM_RDR_rr_1sd.pdf", width = 12, height = 12)
ggplot(summary_F, aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2), color = "black")+ 
  geom_point(size = 8) + 
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], name = "Effect size", ) +
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  #ylim(-0.27, 0.369) +
  labs( x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_markdown(size = 20), legend.position = c(0.2, 0.87), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

# Graph relative rate of model enrichment for all vars - mean and 2*sd
pdf("output/figures/glms/GLM_RDR_rr_2sd.pdf", width = 12, height = 12)
ggplot(summary_F, aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-2*sd_rr_log2, ymax=mean_rr_log2+2*sd_rr_log2), color = "black") +
  geom_point(size = 8) +  
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], name = "Effect size", ) +
  facet_grid(Variable ~ ., scales = "free_y", space = "free_y") +
  #ylim(-0.27, 0.369) +
  labs( x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + theme(axis.text.y = element_markdown(size = 20), legend.position = c(0.2, 0.87), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()


# Range of RDR values
RDR_limits <- range(summary_F$median_pos_RDR, na.rm = TRUE)

# Graph relative rate of model enrichment for all vars - mean and sd
pdf("output/figures/glms/GLM_RDR_rr_abiotic.pdf", width = 8.8398, height = 5.75)
ggplot(summary_F[which(summary_F$Variable == "Abiotic"),], aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2), color = "black") +
  geom_point(size = 8) + 
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], limits = RDR_limits, name = "Effect size", breaks = c(0.04, 0.10, 0.16)) +
  ylim(-0.38, 0.49) +
  labs(title="Abiotic Variables", x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + 
  theme(axis.text.y = element_markdown(size = 22), plot.title = element_text(hjust = 0.5, margin = margin(b = -1)), 
  legend.title = element_text(size = 24), legend.text = element_text(size = 20), legend.position = c(0.22, 0.68), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
pdf("output/figures/glms/GLM_RDR_rr_abiotic_nolegend.pdf", width = 8.8401, height = 5.75)
ggplot(summary_F[which(summary_F$Variable == "Abiotic"),], aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2), color = "black") +
  geom_point(size = 8) + 
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], limits = RDR_limits, name = "Effect size", breaks = c(0.04, 0.10, 0.16)) +
  ylim(-0.38, 0.49) +
  labs(title="Abiotic Variables", x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + 
  theme(axis.text.y = element_markdown(size = 22), plot.title = element_text(hjust = 0.5, margin = margin(b = -1)), legend.position = "none")
dev.off()

# Graph relative rate of model enrichment for all vars - mean and sd
pdf("output/figures/glms/GLM_RDR_rr_biotic.pdf", width = 10, height = 9.5)
ggplot(summary_F[which(summary_F$Variable == "Biotic"),], aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2), color = "black") +
  geom_point(size = 8) + 
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], limits = RDR_limits, name = "Effect size", breaks = c(0.04, 0.10, 0.16)) +
  ylim(-0.38, 0.49) +
  labs(title="Biotic Variables", x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + 
  theme(axis.text.y = element_markdown(size = 22), plot.title = element_text(hjust = 0.5, margin = margin(b = -1)), 
  legend.title = element_text(size = 28), legend.text = element_text(size = 24), legend.position = c(0.77, 0.17), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()
pdf("output/figures/glms/GLM_RDR_rr_biotic_nolegend.pdf", width = 10, height = 9.5)
ggplot(summary_F[which(summary_F$Variable == "Biotic"),], aes(x = reorder(variable_full_name_italics, mean_rr_log2), y = mean_rr_log2, color = median_pos_RDR)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
  geom_errorbar(aes(ymin=mean_rr_log2-sd_rr_log2, ymax=mean_rr_log2+sd_rr_log2), color = "black") +
  geom_point(size = 8) + 
  scale_color_gradientn(colours = brewer.pal(9, "Blues")[2:9], limits = RDR_limits, name = "Effect size", breaks = c(0.04, 0.10, 0.16)) +
  ylim(-0.38, 0.49) +
  labs(title="Biotic Variables", x = "", y = expression(Log[2]*(RR~Enrichment))) +
  coord_flip() +
  theme_linedraw(base_size=30) + 
  theme(axis.text.y = element_markdown(size = 22), plot.title = element_text(hjust = 0.5, margin = margin(b = -1)), legend.position = "none")
dev.off()

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Subset LFMM
summary_sigSNPs_sub <- summary_sigSNPs[,c(1,14,15)]

# Graph GLM RDR vs LFMM z-score metric
effect_GLM_LFMM <- left_join(effectsize_sum, summary_sigSNPs_sub)
effect_GLM_LFMM$top <- c(rep("C", 7), "T", rep("C", 16), "T", "C", "C")

# Graph relative rate of model enrichment for all vars - mean and 2*sd
pdf("output/figures/glms/GLM_RDR_LFMM.pdf", width = 11, height = 8)
ggplot(effect_GLM_LFMM, aes(x = median_pos_RDR, y = abs_zscore_quant_99.9, color = Variable, shape = top)) + 
  geom_point(size = 8) + 
  theme_linedraw(base_size=30) 
dev.off()

