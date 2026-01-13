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
o <- foreach(nL = unique(endemism.merge.filt.SNS$nLocales_poly), .combine="rbind", .errorhandling = "remove"){

    prop <- sum(endemism.merge.filt.SNS$col%like%"missense") / dim(endemism.merge.filt.SNS)[1]
    prop_se <- sqrt(prop * (1-prop)/dim(endemism.merge.filt.SNS)[1])
}

    

    data.table(nLocales_poly=nL,
                mod=c("propNS", "dist", "maf", "score"),
                nSNPs=dim(endemism.merge.filt.SNS)[1],
                p=c(NA, unlist(lapply(mods, function(x) x$p.value))),
                lci=c(prop-1.96*prop_se, unlist(lapply(mods, function(x) x$conf.int[1]))),
                uci=c(prop+1.96*prop_se, unlist(lapply(mods, function(x) x$conf.int[2]))),
                df=c(NA, unlist(lapply(mods, function(x) x$parameter))),
                t=c(NA, unlist(lapply(mods, function(x) x$statistic))),
                delta=c(prop, unlist(lapply(mods, function(x) x$estimate[1]-x$estimate[2]))))
                


# ================================================================================== #

# Graph
pdf("output/figures/endemism/Ncan_endemism.pdf", width = 8, height = 8)
ggplot(endemism.merge.filt.SNS, aes(x=nLocales_poly, )) + 
  geom_line() + 
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()