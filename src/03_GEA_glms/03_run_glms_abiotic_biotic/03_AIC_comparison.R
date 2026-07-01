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
out_dir <- paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
file_names = as.list(dir(path = 'data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/', pattern = "glm.model.collated.*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
glm.model.collated =  foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(glm.model.collated)

# ================================================================================== #

# Save merged data
save(glm.model.collated, file = paste("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.Rdata"))

# ================================================================================== #
# ================================================================================== #

# Test methods with just real
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.model.collated.real.Rdata")

# Identify which column is the minimum
glm.model.collated.real$minAIC <- names(glm.model.collated.real[, 8:13])[apply(glm.model.collated.real[, 8:13], MARGIN = 1, FUN = which.min)]
glm.model.collated.real$minAIC_value <- do.call(pmin, c(glm.model.collated.real[, 8:13], na.rm = TRUE))

# Calc delta AIC between AIC_dem_both_int and model with min AIC
glm.model.collated.real <- glm.model.collated.real %>% mutate(deltaAIC = AIC_dem_both_int-minAIC_value)

# Summarize
glm.real.sum <- glm.model.collated.real %>% group_by(minAIC) |> tally()
glm.real.sum <- glm.real.sum %>% mutate(perc = glm.real.sum$n/dim(glm.model.collated.real)[1]*100)

# Graph delta AIC
pdf("output/figures/GEA/glms/Delta_AIC.pdf", width = 14, height = 6)
ggplot(glm.model.collated.real, aes(x = deltaAIC)) + 
 geom_density(alpha = 0.7, lwd = 0.5) + 
  theme_bw(base_size=30)
dev.off()

pdf("output/figures/GEA/glms/Int_pval.pdf", width = 14, height = 6)
ggplot(glm.model.collated.real, aes(x = p_int)) + 
 geom_density(alpha = 0.7, lwd = 0.5) + 
  theme_bw(base_size=30)
dev.off()
