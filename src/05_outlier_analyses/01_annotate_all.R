# Annotate all SNPs

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
#install.packages(c('data.table', 'tidyverse', 'dplyr', 'geodist', 'foreach', 'doMC'))
library(data.table)
library(tidyverse)
library(dplyr)
library(geodist)
library(foreach)
library(doMC)

# Load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# Register
registerDoMC(20)

# ================================================================================== #

# Specify arguments
args = commandArgs(trailingOnly=TRUE)
i <- as.numeric(args[1])

message("Partition: ", i)

# ================================================================================== #

# Generate output directories

# Data directory
out_dir <- paste("data/processed/outlier_analyses/annotate_all")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Metadata
meta <- read.csv("guide_files/Populations_metadata.csv")
names(meta)[1] <- "sample"

# ================================================================================== #

# Format data

# Extract SNP data from GDS - 14,897,468 SNPs
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Set filter - 14,897,468 SNPs
seqSetFilter(genofile, snp.dt$variant.id)

# ================================================================================== #

# Split snp.dt into 500 chunks
# i.e., create a bin column in snp.dt that goes from 1 to 500 which can be used to subset snp.dt

numJobs <- 500

# Create bins
snp.dt[,bin:=rep(1:numJobs, each=ceiling(dim(snp.dt)[1]/numJobs))[1:dim(snp.dt)[1]]]

# ================================================================================== #

# Iterate through chunks and get allele frequencies

# Check bin
head(snp.dt[bin==i])

# Filter for a specific bin
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=snp.dt[bin==i]$variant.id, sample.id=meta$Site.Code)

# ================================================================================== #

# Get annotations
tmp <- seqGetData(genofile, "annotation/info/ANN")
# Num annotations per variant
len1 <- tmp$length
# Annotation data
len2 <- tmp$data

# Create df of annotations and variant ids
snp.dt1 <- data.table(
            len=rep(len1, times=len1),
            ann=len2,
            id=rep(snp.dt[bin==i]$variant.id, times=len1))

# Extract data between the 2nd and third | symbol
snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

# Collapse additional annotations to original SNP vector length
snp.dt1.an <- snp.dt1[,list(n_ann=length(class), col= paste(class, collapse=","), gene=paste((gene), collapse=",")),
                       list(variant.id=id)]

# Split the col and gene columns if there is a comma and only keep the first; also replace empty with NA
snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

# Merge
o <- merge(snp.dt[bin==i], snp.dt1.an, by="variant.id")

# Save output
save(o, file= paste("data/processed/outlier_analyses/annotate_all/chunk_", i, ".Rdata", sep=""))