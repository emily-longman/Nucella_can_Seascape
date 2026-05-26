# Analyze windows

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'doMC'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/outlier_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load data - outliers for pH and Mcali based on beating POD BF threshold
bf.ph.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header = T)
bf.McaliIntThk.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", header = T)

# Load annotation
load("data/processed/outlier_analyses/Ncan.pooldata.annotations.RData")

# ================================================================================== #

# Summarize each list based on annotation
ph_ann_sum <- bf.ph.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(ph = n)
McaliIntThk_ann_sum <- bf.McaliIntThk.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(McaliIntThk = n)

# Join data
ann <- full_join(ph_ann_sum, McaliIntThk_ann_sum)

# Summarize full SNP list
all_ann <- o.merge.filt %>% count(col) %>% rename(Annotation = col, genome = n)

# Join data
ann <- full_join(ann, all_ann)
ann[is.na(ann)] <- 0

#--------------------------------------------------------------------------------

# Subset data for common var
ann.focal <- c("missense_variant", "synonymous_variant", "intron_variant", "intergenic_region", "upstream_gene_variant", "downstream_gene_variant")
ann.sub <- ann %>% filter(Annotation %in% ann.focal)

# Loop through annotations and perform Fishers exact test

ftests <- foreach(i=ann.focal, .combine = "rbind", .errorhandling = "remove")%do%{
  
  # Create matrix for significance of focal annotation
  ann.tmp <- matrix(c(
            ann$ph[which(ann$Annotation==i)],
            ann$McaliIntThk[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$ph[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$McaliIntThk[which(ann$Annotation==i)]),
            nrow = 2)

  # Fishers exact test
  ftest_i <- fisher.test(ann.tmp)

  # Create data frame with output
  data.frame(
      class = i,
      p.fet = ftest_i$p.value,
      OR = ftest_i$estimate,
      lci = ftest_i$conf.int[1],
      uci = ftest_i$conf.int[2]
    )
}

#--------------------------------------------------------------------------------

# Change class to factor
ftests$class <- factor(ftests$class, levels = ann.focal)

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test.pdf", width = 6, height = 6)
ggplot(ftests, aes(y=log2(OR), x = class)) + 
  geom_point(size=3.5) + 
  geom_linerange(aes(ymin = log2(lci), ymax = log2(uci))) + 
  geom_hline(yintercept = 0, col="black", linetype="dashed") + xlab("") +
  theme_bw(base_size=16) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()
