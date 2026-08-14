# Graph endemism analysis

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
#install.packages(c('data.table', 'tidyverse', 'dplyr', 'foreach', 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(dplyr)
library(foreach)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/endemism")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load filtered data
load("data/processed/endemism/endemism.merge.filt.RData")

# ================================================================================== #

# Filter for only S and NS
endemism.merge.filt.SNS <- endemism.merge.filt %>% filter(col == "synonymous_variant" | col == "missense_variant")

# Summarize
o <- foreach(nL=unique(endemism.merge.filt.SNS$nLocales_poly), .combine="rbind", .errorhandling="remove")%do%{

    tmp <- endemism.merge.filt.SNS %>% filter(nLocales_poly == nL)

    propNS <- sum(tmp$col%like%"missense_variant") / dim(endemism.merge.filt.SNS)[1]
    propNS_se <- sqrt(propNS * (1-propNS)/dim(endemism.merge.filt.SNS)[1])

    data.table(
       nLocales_poly = nL,
       nSNP=dim(tmp)[1],
       propNS = propNS,
       propNS_se = propNS_se
    )
}

# ================================================================================== #

# Graph
pdf("output/figures/endemism/Ncan_endemism.pdf", width = 8, height = 8)
ggplot(o, aes(x=nLocales_poly, propNS)) + 
  geom_line() + 
  ylab("Proportion of\nmissense coding SNPs") +
  geom_hline(aes(yintercept=0), linetype="dashed") +
  theme_classic(base_size = 20)
dev.off()