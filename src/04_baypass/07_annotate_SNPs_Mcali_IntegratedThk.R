# This script will annotate SNPs of interest.

# ================================================================================== #

# Set path as main Github repo
# Install and load package
install.packages(c('rprojroot'))
library(rprojroot)
# Specify root path
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

#--------------------------------------------------------------------------------

# Load packages
install.packages(c('data.table', 'tidyverse', 'foreach'))
library(data.table)
library(tidyverse)
library(foreach)

# Load SeqArray
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install("SeqArray")
library(SeqArray)

#--------------------------------------------------------------------------------

# Convert VCF to GDS

# Load the VCF file
#vcf.fn <- "data/processed/outlier_analyses/snpeff/N.canaliculata_pops_SNPs_annotate.vcf"
# Parse the header
#seqVCF_Header(vcf.fn)
# Convert VCF to GDS
#seqVCF2GDS(vcf.fn, storage.option="ZIP_RA", "data/processed/outlier_analyses/snpeff/N.canaliculata_annotated_SNPs.gds")

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

#--------------------------------------------------------------------------------

# Load SNPs of interest
bf.McaliIntThk.mean.sum.outliers <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.csv", header=T)

#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest

annotation <- foreach(i=1:dim(bf.McaliIntThk.mean.sum.outliers)[1], .combine = "rbind", .errorhandling = "remove")%do%{

  message(i)
  # Reset filter
  seqResetFilter(genofile)
  # Extract SNP_id for SNP i
  tmp.i = bf.McaliIntThk.mean.sum.outliers[i,]$SNP_id
  # Extract snp.dt information for SNP i
  pos.tmp = snp.dt %>% filter(SNP_id %in% tmp.i) %>% .$id
  # Set filter for SNP i
  seqSetFilter(genofile, variant.id = pos.tmp)
  # Extract annotation
  ann_data <- seqGetData(genofile, "annotation/info/ANN")$data
  # Identify if multiple annotation
  L = length(ann_data)

  # Loop through annotations for SNP i
  annotate.list =
  
  foreach(k=1:L, .combine = "rbind")%do%{

    tmp = ann_data[k] 
    tmp2= str_split(tmp, "\\|")
  
    data.frame(
      id=pos.tmp,
      SNP_id = tmp.i,
      annotation.id=k,
      Allele = tmp2[[1]][1],
      Annotation = tmp2[[1]][2],
      Annotation_Impact = tmp2[[1]][3],
      Gene_Name = tmp2[[1]][4],
      Gene_ID = tmp2[[1]][5],
      Feature_Type = tmp2[[1]][6],
      Feature_ID = tmp2[[1]][7],
      Transcript_BioType = tmp2[[1]][8],
      Rank = tmp2[[1]][9],
      HGVS.c = tmp2[[1]][10],
      HGVS.p = tmp2[[1]][11],
      cDNA.pos.cDNA.length = tmp2[[1]][12],
      CDS.pos.CDS.length = tmp2[[1]][13],
      AA.pos.AA.length = tmp2[[1]][14],
      Distance = tmp2[[1]][15]
      )
  }
return(annotate.list)
}

#--------------------------------------------------------------------------------

# Join annotation and SNP information
bf.McaliIntThk.mean.sum.outliers.annotated <- left_join(bf.McaliIntThk.mean.sum.outliers, annotation, by = join_by(SNP_id))

# Write output
write.csv(bf.McaliIntThk.mean.sum.outliers.annotated, "data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", row.names = F, quote = F)
