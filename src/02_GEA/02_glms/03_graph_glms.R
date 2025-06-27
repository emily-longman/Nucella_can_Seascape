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

# Generate output directory

out_dir <- paste("data/processed/GEA/glms/glms_output")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

out_fig_dir <- paste("output/figures/GEA/glms")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_window_analysis/', pattern = "GLM_100perm_Bio-Oracle_chunk_*"))
file_names = lapply(file_names, function(x) paste0('data/processed/GEA/glms/glms_window_analysis/', x))

# Read all the files and add a column with the chunk
glm.model.collated =  
foreach(i=file_names, .combine="rbind")%do%{  
message(i) 
chunk = i
o = get(load(i))
o %>% mutate(chunk = file_names) %>% mutate(chunk = str_remove(chunk, pattern = "data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_"))
} 

# Save merged data
save(glm.model.collated, file = "data/processed/GEA/glms/glms_output/glm.model.collated.Rdata")

# ================================================================================== #

# Add SNP_id column
glm.model.collated %>% mutate(SNP_id = paste(chr, pos, sep = "_"))


# ================================================================================== #

# Graph results
pdf("output/figures/GEA/glms/glm_perm.pdf", width = 8, height = 8)
ggplot(data=glm.model.collated, aes(x=variable, y=-log10(p_lrt), group=perm, color = perm)) + 
geom_boxplot()
dev.off()

pdf("output/figures/GEA/glms/glm_perm.pdf", width = 8, height = 8)
ggplot(data=glm.model.collated[glm.model.collated$data == "real"], aes(x=-log10(p_lrt), y=N, group=perm, color=as.factor(perm!=0))) + 
geom_line()
dev.off()



# not finished - took from JCBN code

pdf("output/figures/morphology/GLM.pdf", width = 6, height = 6)
GLM_test %>%
  group_by(perm==0, variant.id) %>%
  summarise( uci = quantile(p_lrt, 0.1, na.rm = T)) %>%
  separate(variant.id, remove = F, into = c("chr_Or", "chr_id", "pos"), sep = "_" ) %>%
  ggplot(aes(
    x=as.numeric(pos),
    y=-log10(uci),
    color=`perm == 0`
  )) + geom_line()
dev.off()

GLM_test %>%
  group_by(perm==0, variant.id) %>%
  summarise( uci = quantile(p_lrt, 0.1, na.rm = T)) %>%
  separate(variant.id, remove = F, into = c("chr_Or", "chr_id", "pos"), sep = "_" ) %>%
  dcast(variant.id+chr_Or+chr_id+pos~`perm == 0`) %>%
  mutate(test = `TRUE`<`FALSE`) %>%
  filter(test == "TRUE")


