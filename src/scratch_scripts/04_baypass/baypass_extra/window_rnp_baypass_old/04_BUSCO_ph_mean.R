# Extract just BUSCO

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
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #


# Read outlier SNPs
outlier.win.SNPs <- read.csv("data/processed/baypass/window_summary/window_100kb_ph_mean_outlier_SNPs.csv", header=T)

# ================================================================================== #

# Load BUSCO table
BUSCO <- fread("data/processed/BUSCO/N_canaliculata/run_mollusca_odb12/full_table.tsv", skip = 2, fill=TRUE)
colnames(BUSCO) <- c("Busco_id", "Status", "Sequence", "Gene_Start", "Gene_End", "Strand", "Score", "Length", "OrthoDB_url", "Description")

# Filter for just complete
BUSCO.complete <- BUSCO %>% filter(Status == "Complete")

# Get SNPs in each complete BUSCO
BUSCO.complete.SNPs <- foreach(busco.i=unique(BUSCO.complete$Busco_id), .combine="rbind", .errorhandling="remove")%do%{
    
    # Extract window
    busco.tmp <- BUSCO.complete[which(BUSCO.complete$Busco_id==busco.i),]

    # Extract SNPs in window from baypass output
    snpdet %>% filter(
        snpdet$chr == busco.tmp$Sequence & 
        snpdet$pos > busco.tmp$Gene_Start &
        snpdet$pos < busco.tmp$Gene_End)
}

# Create SNP_id columns
BUSCO.complete.SNPs <- BUSCO.complete.SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
outlier.win.SNPs <- outlier.win.SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# Find overlap
overlap_SNP_id <- intersect(BUSCO.complete.SNPs$SNP_id, outlier.win.SNPs$SNP_id)

# Extract SNPs in outlier windows in complete BUSCOs (2,896 SNPs)
outlier.win.SNPs.busco <- outlier.win.SNPs %>% filter(SNP_id %in% overlap_SNP_id)

# ================================================================================== #

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_annotated_SNPs.gds")

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        id=seqGetData(genofile, "variant.id")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

# Extract annotation data for each SNP of interest

annotation <- foreach(i=1:dim(outlier.win.SNPs.busco)[1], .combine = "rbind", .errorhandling = "remove")%do%{

  message(i)
  # Reset filter
  seqResetFilter(genofile)
  # Extract SNP_id for SNP i
  tmp.i = outlier.win.SNPs.busco[i,]$SNP_id
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
outlier.win.SNPs.busco.annotated <- left_join(outlier.win.SNPs.busco, annotation, by = join_by(SNP_id), relationship = "many-to-many")

# Write output
write.csv(outlier.win.SNPs.busco.annotated, "data/processed/baypass/window_summary/outlier.win.SNPs.busco.annotated.ph_mean.csv", row.names = F, quote = F)

