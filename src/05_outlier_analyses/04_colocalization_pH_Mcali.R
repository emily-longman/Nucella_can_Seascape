# Co-localization of BF of GWAS and GEA

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
install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'groupdata2', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(groupdata2)
library(poolfstat)
library(RColorBrewer)
library(viridis)
library(SeqArray)

# ================================================================================== #

# Load data

# Load GEA data
# pH Data
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
# M. cali thk
load("data/processed/baypass/biotic/bf.Mcali.IntThk.sum.Rdata")

# Load thresholds
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")
bf.POD.thr.ph <- bf.POD.thr
load("data/processed/baypass/biotic/Mcali_IntegratedThk_POD_thr.Rdata")
bf.POD.thr.McaliThk <- bf.POD.thr

# ================================================================================== #

# Compare pH and Mcali

# Rename cols
bf.ph.mean.sum <- bf.ph.mean.sum %>% rename(bf_db.mean.ph = bf_db.mean, bf_db.median.ph = bf_db.median, bf_db.var.ph = bf_db.var, 
eBPis.mean.ph = eBPis.mean, eBPis.median.ph = eBPis.median, eBPis.var.ph = eBPis.var)

# pH and Mcali
bf.ph.Mcali.sum <- left_join(bf.ph.mean.sum, bf.McaliIntThk.mean.sum, by=c("chr", "pos", "allele1", "allele2", "MRK"))

# Graph ph vs Mcali
pdf("output/figures/outlier_analyses/BF_ph_vs_Mcali.pdf", width = 8, height = 8)
ggplot(bf.ph.Mcali.sum, aes(x=bf_db.mean.ph, y=bf_db.mean)) +
  labs(x = "BF mean pH", y = "BF M. californianus thickness") +
  geom_point(alpha=0.6) +
  geom_vline(xintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  theme_classic(base_size = 30)
dev.off()

pdf("output/figures/outlier_analyses/BF_ph_vs_Mcali_posBF.pdf", width = 8.5, height = 8.14)
ggplot(bf.ph.Mcali.sum[which(bf.ph.Mcali.sum$bf_db.mean.ph>0 & bf.ph.Mcali.sum$bf_db.mean>0),], aes(x=bf_db.mean.ph, y=bf_db.mean)) +
  labs(x = "BF mean pH", y = expression(paste("BF ", italic("M. californianus"), " thickness"))) +
  geom_point(alpha=0.6, col="#4a4949") +
  geom_vline(xintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  geom_point(data = bf.ph.Mcali.sum[which(bf.ph.Mcali.sum$bf_db.mean.ph>14 & bf.ph.Mcali.sum$bf_db.mean>12),], aes(x=bf_db.mean.ph, y=bf_db.mean), col="black",) +
  theme_linedraw(base_size = 30)
dev.off()
pdf("output/figures/outlier_analyses/BF_ph_vs_Mcali_posBF2.pdf", width = 7.8, height = 7.8)
ggplot(bf.ph.Mcali.sum[which(bf.ph.Mcali.sum$bf_db.mean.ph>0 & bf.ph.Mcali.sum$bf_db.mean>0),], aes(x=bf_db.mean.ph, y=bf_db.mean)) +
  labs(x = "BF mean pH", y = "BF mussel shell thickness") +
  geom_point(alpha=0.6, col="#4a4949") +
  geom_vline(xintercept=bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)], col="red") +
  geom_hline(yintercept=bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)], col="red") +
  geom_point(data = bf.ph.Mcali.sum[which(bf.ph.Mcali.sum$bf_db.mean.ph>14 & bf.ph.Mcali.sum$bf_db.mean>12),], aes(x=bf_db.mean.ph, y=bf_db.mean), col="black",) +
  theme_linedraw(base_size = 30)
dev.off()

# ================================================================================== #

# Which SNPs are closest to cross
topSNP <- bf.ph.Mcali.sum[
    which(bf.ph.Mcali.sum$bf_db.mean.ph > bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)]-1.0 
          & bf.ph.Mcali.sum$bf_db.mean > bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)]-1.0),]

# Extend range a bit (BF ph > 9.37 and BF Mcali > 8.5)
topSNP_extra <- bf.ph.Mcali.sum[
    which(bf.ph.Mcali.sum$bf_db.mean.ph > bf.POD.thr.ph$bf_db.mean[which(bf.POD.thr.ph$thr==0.999)]-5.0 
          & bf.ph.Mcali.sum$bf_db.mean > bf.POD.thr.McaliThk$bf_db.mean[which(bf.POD.thr.McaliThk$thr==0.999)]-5.0),]

# ================================================================================== #

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

# Make SNP_id column for outliers
topSNP_extra <- topSNP_extra %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest

annotation <- foreach(i=1:dim(topSNP_extra)[1], .combine = "rbind", .errorhandling = "remove")%do%{

  message(i)
  # Reset filter
  seqResetFilter(genofile)
  # Extract SNP_id for SNP i
  tmp.i = topSNP_extra[i,]$SNP_id
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
topSNP_extra.annotated <- left_join(topSNP_extra, annotation, by = join_by(SNP_id))

# Write output
write.csv(topSNP_extra.annotated, "data/processed/outlier_analyses/topSNP_ph_Mcali_extra.annotated.csv", row.names = F, quote = F)
