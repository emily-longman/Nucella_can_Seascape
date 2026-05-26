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

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load data - outliers for pH and Mcali based on beating POD BF threshold
bf.ph.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header = T)
bf.McaliIntThk.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", header = T)

# ================================================================================== #

# Summarize each list based on annotation
ph_ann_sum <- bf.ph.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(ph = n)
McaliIntThk_ann_sum <- bf.McaliIntThk.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(McaliIntThk = n)

# Join data
ann <- full_join(ph_ann_sum, McaliIntThk_ann_sum)
ann[is.na(ann)] <- 0
row.names(ann) <- ann$Annotation

#--------------------------------------------------------------------------------

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_annotated_SNPs.gds")

#--------------------------------------------------------------------------------

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        id=seqGetData(genofile, "variant.id")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Summarize specific fields

seqResetFilter(genofile)
# Extract SNP_id for SNP i
tmp.i = snp.dt[1,]$SNP_id
pos.tmp = snp.dt %>% filter(SNP_id %in% tmp.i) %>% .$id

seqSetFilter(genofile, variant.id = pos.tmp)


ann_data <- seqGetData(genofile, "annotation/info/ANN")$data



  foreach(k=1:L, .combine = "rbind")%do%{
    tmp = ann_data[k] 
    tmp2= str_split(tmp, "\\|")
  
    data.frame(
      Annotation = tmp2[[1]][2],
      )
      }

### FIRST PASS

# Subset for just common ann
ann.sub <- ann[c("3_prime_UTR_variant", "5_prime_UTR_variant", "downstream_gene_variant", "intergenic_region", "intron_variant", "missense_variant", "synonymous_variant", "upstream_gene_variant"),]

pdf("output/figures/baypass/window_summary/mosaic_plot.pdf", width = 12, height = 6)
mosaicplot(ann.sub, color = TRUE)
dev.off()

# Fishers exact test (Monte Carlo simulation with 2,000 simulations)
ftest <- fisher.test(ann.sub, simulate.p.value = TRUE, B = 2000)

