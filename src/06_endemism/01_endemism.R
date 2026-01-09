# Endemism analysis
# Modified from L. Proud; Original from A. Bergland

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
#install.packages(c('SeqArray', 'data.table', 'tidyverse', 'dplyr', 'geodist', 'foreach', 'doMC'))
library(SeqArray)
library(data.table)
library(tidyverse)
library(dplyr)
library(geodist)
library(foreach)
library(doMC)

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
out_dir <- paste("data/processed/endemism")
if (!dir.exists(out_dir)) {dir.create(out_dir)}

# ================================================================================== #

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs.annotate.gds")

# Metadata
meta <- read.csv("guide_files/Populations_metadata.csv")

# ================================================================================== #

# Format data

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        variant.id=seqGetData(genofile, "variant.id"),
        allele=seqGetData(genofile, "allele")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))


#head(snp.dt[,list(minPos=min(pos), maxPos=max(pos)), list(chr)])

# ================================================================================== #

# Subset to sites with only 2 alleles
snp.dt <- snp.dt[nAlleles==2]
seqSetFilter(genofile, snp.dt$variant.id) # 14,691,690 SNPs

# ================================================================================== #

# Chunk snp.dt - i.e., create a bin column that goes from 1 to 500 which can be used to subset snp.dt

numJobs <- 500

# Create bins
snp.dt[,bin:=rep(1:numJobs, each=ceiling(dim(snp.dt)[1]/numJobs))[1:dim(snp.dt)[1]]]


####
# Load pooldata object
load("data/raw/pooldata/pooldata.RData")

# Extract SNP info for all SNPs
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) ->
  pooldata.snp.info

# Rename columns
names(pooldata.snp.info)[1:2] = c("chr","pos")

# Make snp_id column
pooldata.snp.info <- pooldata.snp.info %>%
  mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Iterate through chunks and get allele frequencies

# Check bin
head(snp.dt[bin==i])

# Filter for a specific bin
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=snp.dt[bin==i]$variant.id, sample.id=meta$Site.Code)

# Get allele frequencies (take every other column so just get minor allele) and depth
ad <- seqGetData(genofile, "annotation/format/AD") %>% .$data
ad.minor <- ad[, seq(2, ncol(ad), by = 2)]
dp <- seqGetData(genofile, "annotation/format/DP")

# Get sample ids and variants
samples <- seqGetData(genofile, "sample.id")
variants <- seqGetData(genofile, "variant.id")

# Confirm lengths match the dimensions
stopifnot(length(samples) == dim(ad.minor)[1])
stopifnot(length(variants) == dim(ad.minor)[2])

# Assign dimnames as a list of two vectors (variant and sample)
dimnames(ad.minor) <- list(sample = samples, variant = variants)
dimnames(dp) <- list(sample = samples, variant = variants)

# Reformat data
tmp.ad.minor <- as.data.table(reshape2::melt(ad.minor))
tmp.dp <- as.data.table(reshape2::melt(dp))

# Setkey
setkey(tmp.ad.minor, sample, variant)
setkey(tmp.dp, sample, variant)

# Merge ad and dp
m <- merge(tmp.ad.minor, tmp.dp)
# Rename columns
setnames(m, c("value.x", "value.y"), c("ad", "dp"))

# Check structure
head(m)

# Merge table with meta
m <- data.table(merge(meta, m, by.x="Site.Code", by.y="sample"))

# Calculate allele frequency
m$freq <- m$ad/m$dp

# ================================================================================== #

# Summarize
setkey(m, Site.Code)

m.ag <- m[,list(nSamps_poly=sum(freq>0 & freq<1, na.rm=T), 
                nSamps_fixed=sum(freq==0 | freq==1, na.rm=T),
                nSamps_missing=sum(is.na(freq)),

                nLocales_poly=length(unique(na.omit(Site.Code[freq>0 & freq<1]))),
                #nSamps2_poly=length(unique(na.omit(sampleId[freq>0 & freq<1]))), #Unlike the flies, I don't have multiple samples per site

                nLocales_fixed=length(unique(na.omit(Site.Code[freq==0 | freq==1]))),
                #nSamps2_fixed=length(unique(na.omit(sampleId[freq==0 & freq==1]))),

                global_af=sum(ad, na.rm=T)/sum(dp, na.rm=T),
                poly_af=mean(freq[freq>0 & freq<1 & !is.na(freq)], na.rm=T),
                poly_ad_mean=mean(ad[freq>0 & freq<1 & !is.na(freq)], na.rm=T)+.00001,
                poly_ad_median=median(ad[freq>0 & freq<1 & !is.na(freq)], na.rm=T)+.00001,
                poly_min_ad=min(ad[freq>0 & freq<1 & !is.na(freq)], na.rm=T)+.00001,
                poly_max_ad=max(ad[freq>0 & freq<1 & !is.na(freq)], na.rm=T)+.00001,
                poly_samps=paste(Site.Code[freq>0 & freq<1 & !is.na(freq)], collapse=";")
                #poly_locales=paste(sort(unique(na.omit(locality[freq>0 & freq<1]))), collapse=";")

                  ),

                  list(variant)]

m.ag <- m.ag[nSamps_poly!=0]

# Check structure
print(str(m.ag))          # Make sure it's a data.table and has expected columns
print(head(m.ag))         # Check first few rows
print(length(m.ag$variant))  # Confirm length matches your expectation

# ================================================================================== #

# NOTE: somehow the numbers get off here between df

# Get annotations
message("Annotations")
seqResetFilter(genofile)
seqSetFilter(genofile, variant.id=m.ag[nSamps_poly!=0]$variant)
s <- seqGetData(genofile, "annotation/info/ANN")
len1 <- tmp$length
len2 <- tmp$data

snp.dt1 <- data.table(
            len=rep(len1, times=len1),
            ann=len2,
            id=rep(m.ag$variant, times=len1))

# Extract data between the 2nd and third | symbol
snp.dt1[,class:=tstrsplit(snp.dt1$ann,"\\|")[[2]]]
snp.dt1[,gene:=tstrsplit(snp.dt1$ann,"\\|")[[4]]]

# Collapse additional annotations to original SNP vector length
snp.dt1.an <- snp.dt1[,list(n=length(class), col= paste(class, collapse=","), gene=paste((gene), collapse=",")),
                       list(variant.id=id)]

    
snp.dt1.an[,col:=tstrsplit(snp.dt1.an$col,"\\,")[[1]]]
snp.dt1.an[,gene:=tstrsplit(snp.dt1.an$gene,"\\,")[[1]]]

m.ag <- merge(m.ag, snp.dt1.an, by.x="variant", by.y="variant.id")

#table(m.ag$col=="stop_gained", m.ag$nSamps_poly)
#m.ag[col=="stop_gained" & nSamps_poly==5]
o <- merge(snp.dt[bin==i], m.ag, by.x="variant.id", by.y="variant")

      
#o <- rbindlist(o)

#system("mkdir /scratch/aob2x/private_v3_slices")

x <- o[nSamps_poly>0]
x

save(x, file= paste("data/processed/endemism/slice", job, ".Rdata", sep=""))
  
##  save(o, file=paste("/scratch/aob2x/private_v5_slices/slice", job, ".Rdata", sep=""))