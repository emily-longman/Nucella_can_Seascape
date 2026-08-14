# Extract xtx outliers

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
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)

# ================================================================================== #

# Load data

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load and merge xtx data
NC.xtx <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/xtx/NC_run", i, "_summary_pi_xtx.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}
colnames(NC.xtx) <- c("chr", "pos", "allele1", "allele2", "MRK", "M_P", "SD_P", "M_XtX", "SD_XtX", "XtXst", "log10.1.pval.", "rep")
# Average xtx across replicate runs
NC.xtx.sum <- NC.xtx %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(XtXst.mean = mean(XtXst), XtXst.median = median(XtXst), log10.1.pval.mean = mean(log10.1.pval.), log10.1.pval.median = median(log10.1.pval.))

# ================================================================================== #

# Identify SNPs undergoing positive selection (i.e., significant SNPs with XtXst higher than the mean) 


# Calculate mean xtx
mean.xtx <- mean(NC.xtx.sum$XtXst.mean)
# Standardize xtx
NC.xtx.sum <- NC.xtx.sum %>% mutate(XtXst.mean.standardize = XtXst.mean - mean.xtx)
# Identify SNPs with XtXst higher than the mean and are significant
NC.xtx.sum.pos <- NC.xtx.sum %>% filter(XtXst.mean.standardize > 0 & log10.1.pval.mean > -log10(0.001))

# Save Outliers
save(NC.xtx.sum.pos, file="data/processed/genomic_offset/NC.xtx.sum.pos.Rdata")

###

# Group based on pos sel, balancing sel and neutrals
NC.xtx.sum <- NC.xtx.sum %>% mutate(group = case_when(
    XtXst.mean.standardize > 0 & log10.1.pval.mean > -log10(0.001) ~ "PosSel", 
    XtXst.mean.standardize < 0 & log10.1.pval.mean > -log10(0.001) ~ "BalSel",
   TRUE ~ "Neutral"))

# Graph xtx
pdf("output/figures/baypass/xtx_colored.pdf", width = 15, height = 8)
ggplot(NC.xtx.sum, aes(y=XtXst.mean, x=chr, col=group)) +
  labs(x = "Position", y = "XtX") + scale_color_viridis_d() +
  geom_point(alpha=0.6) + theme_classic(base_size = 20) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
dev.off()
