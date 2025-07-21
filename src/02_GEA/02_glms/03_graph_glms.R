# Graph glms output

# Clear memory
rm(list=ls())

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('data.table', 'tidyverse', 'plyr', 'foreach', 'ggplot2'))
library(data.table)
library(tidyverse)
library(plyr)
library(foreach)
library(ggplot2)

# ================================================================================== #

# Generate output directories

out_dir <- paste("data/processed/GEA/glms/glms_output")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

out_fig_dir <- paste("output/figures/GEA/glms")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_test/', pattern = "GLM_100perm_Bio-Oracle_chunk_*"))
file_names = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_chunk_analysis_test/', x))))

# Read all the files and add a column with the chunk
glm.model.collated =  
foreach(i=file_names, .combine="rbind")%do%{  
message(i) 
chunk = i
o = get(load(i))
o %>% mutate(chunk = file_names) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_chunk_analysis_test/GLM_100perm_Bio-Oracle_"))
} 

# Save merged data
save(glm.model.collated, file = "data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")

load("data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")

# ================================================================================== #

# Add SNP_id column
glm.model.collated <- glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Split data by environmental variable

# Get list of environmental var names
unique(glm.model.collated$variable) -> enviro_vars_names

# Group by environmental variable
glm.model.collated.group <- glm.model.collated %>% group_by(variable)

# Split by environmental variable
glm.model.collated.split <- group_split(glm.model.collated.group)

# Extract environmental var names
group_keys(glm.model.collated.group) -> enviro_vars_names

# ================================================================================== #

# Extract real data and permutation data
real_data <- glm.model.collated %>% filter(perm == 0)
perm_data <- glm.model.collated %>% filter(perm > 0)


# How do I summarize the permutation data - mean of p_lrt? quantiles? Then how do I compare this to the real data?

# Not correct below
ratios <- perm_data %>%
left_join(real_data, by = "variable", suffix = c("_perm", "_real")) %>%
mutate(rr_ratio = rr_real / rr_perm) %>%
select(variable, test_code_perm, rr_ratio)

# Summarize
summary <- glm.model.collated %>% group_by(variable) %>% 
summarise(mean = mean(p_lrt),
            sd = sd(p_lrt),
            n = n(),  # Number of observations
            se = sd / sqrt(n),  # Calculate standard error
            ci_low = mean - 1.96 * se,
            ci_high = mean + 1.96 * se,
            quantile_0.05 = quantile(p_lrt, 0.05),
            quantile_0.95 = quantile(p_lrt, 0.95),
            quantile_0.5 = quantile(p_lrt, 0.5))


# ================================================================================== #


# Graph results
#pdf("output/figures/GEA/glms/glm_perm.pdf", width = 8, height = 8)
#ggplot(data=glm.model.collated, aes(x=variable, y=-log10(p_lrt), group=perm, color = perm)) + 
#geom_boxplot()
#dev.off()

#pdf("output/figures/GEA/glms/glm_perm.pdf", width = 8, height = 8)
#ggplot(data=glm.model.collated[glm.model.collated$data == "real"], aes(x=-log10(p_lrt), y=N, group=perm, color=as.factor(perm!=0))) + 
#geom_line()
#dev.off()


pdf("output/figures/GEA/glms/glm_perm.pdf", width = 8, height = 8)
glm.model.collated.split[1] %>%
  group_by(perm == 0, SNP_id) %>%
  summarise(uci = quantile(p_lrt, 0.1, na.rm = T)) %>%
  separate(SNP_id, remove = F, into = c("chr", "pos"), sep = "_" ) %>%
  ggplot(aes(
    x=chr,
    y=-log10(uci),
    color=`perm == 0`
  )) + geom_line()
dev.off()
