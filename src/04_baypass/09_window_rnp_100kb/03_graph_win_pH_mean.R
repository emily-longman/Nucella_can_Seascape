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
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/baypass/window_summary/window_100kb_chunk_analysis_ph_mean/")
file_names = as.list(dir(path = path, pattern = "window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/baypass/window_summary/window_100kb_chunk_analysis_ph_mean/"), x))))

# Check number of files
length(file_names_v)

# Read all the files and add a column with the chunk
win.out =  
foreach(i=file_names_v, .combine="rbind", .errorhandling = "remove")%do%{  
    # State which file loading
    message(i)
    # Load file
    o = get(load(i))
}

# Check structure
str(win.out)

# ================================================================================== #

# Read in SNP data from Baypass
snpdet <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snpdet) <- c("chr", "pos", "allele1", "allele2")
# Make unique list of chr names
snpdet.chr <- unique(snpdet$chr)

# ================================================================================== #

# Make sure windows are ordered in same chr list as snpdet
win.out.order <- win.out[order(factor(win.out$chr, levels = snpdet.chr)),]

# Note: number of chr between win.out.order and snpdet don't match becuase several chr failed the filter when generating the windows

# ================================================================================== #

# Save merged data
save(win.out.order, file = "data/processed/baypass/window_summary/window_analysis_ph_mean.RData")
load("data/processed/baypass/window_summary/window_analysis_ph_mean.RData")

# ================================================================================== #

# Use the POD threshold to come up with p-val

# Load mean bf data from 5 baypass runs
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
# Load POD thresholds
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")

# Create the ECDF (empirical cumulative distribution function) function
my_ecdf <- ecdf(bf.ph.mean.sum$bf_db.mean)

# Find the probability for a given value
probability <- my_ecdf(bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)])

# Calc pr.i as the opposite of the probability
pr.i <- c(1-probability)

# ================================================================================== #

# Create unique Chromosome number
win.out.order.chr.unique <- unique(win.out.order$chr)
win.out.order$chr.unique <- as.numeric(factor(win.out.order$chr, levels = win.out.order.chr.unique))

# Graph rnp p

# Graph rnp geompoint
pdf("output/figures/baypass/window_summary/baypass_window_ph_mean_rnpPOD_geompoint.pdf", width = 12, height = 6)
ggplot(win.out.order, aes(y=-log10(rnp.binom.POD), x=chr.unique)) + 
  geom_point(alpha=0.8, size=1.6) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=26) + theme(legend.position = "none")
dev.off()

# Graph rnp geomline
pdf("output/figures/baypass/window_summary/baypass_window_ph_mean_rnpPOD_geomline_wider.pdf", width = 12, height = 3)
ggplot(win.out.order, aes(y=-log10(rnp.binom.POD), x=chr.unique)) + ylab("100kb Window\nEnrichment") + xlab("Scaffold") +
  geom_line( ) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=18) + theme(legend.position = "none") + theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()

# Graph Windows
pdf("output/figures/baypass/window_summary/baypass_window_ph_mean_rnpPOD_geomline_bywindows.pdf", width = 12, height = 3)
ggplot(win.out.order, aes(y=-log10(rnp.binom.POD), x=win)) + ylab("100kb Window\nEnrichment") + xlab("               Window") +
  geom_line( ) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=18) + theme(legend.position = "none") #+ theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()

# ================================================================================== #

# Extract outliers
win.out.order.outliers <- win.out.order %>% filter(-log10(rnp.binom.POD) > -log10(pr.i))

# Save outliers
write.csv(win.out.order.outliers, "data/processed/baypass/window_summary/window_100kb_analysis_ph_mean_outliers.csv", row.names=FALSE)

# ================================================================================== #

# Identify SNPs in the outlier windows -> 45,095 SNPs
outlier.win.SNPs <- foreach(win.i=unique(win.out.order.outliers$win), .combine="rbind", .errorhandling="remove")%do%{
    
    # Extract window
    win.tmp <- win.out.order.outliers[which(win.out.order.outliers$win==win.i),]

    # Extract SNPs in window from baypass output
    baypass.tmp <- bf.ph.mean.sum %>% filter(
        bf.ph.mean.sum$chr == win.tmp$chr & 
        bf.ph.mean.sum$pos >= win.tmp$pos_min &
        bf.ph.mean.sum$pos <= win.tmp$pos_max)
}

# ================================================================================== #

# Add SNP_id column
outlier.win.SNPs <- outlier.win.SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# Save outlier SNPs
write.csv(outlier.win.SNPs, "data/processed/baypass/window_summary/window_100kb_ph_mean_outlier_SNPs.csv", row.names=FALSE)
outlier.win.SNPs <- read.csv("data/processed/baypass/window_summary/window_100kb_ph_mean_outlier_SNPs.csv", header=T)

#--------------------------------------------------------------------------------

# Open the GDS file
genofile <- seqOpen("data/processed/outlier_analyses/snpeff/N.canaliculata_annotated_SNPs.gds")

# Extract SNP data from GDS
snp.dt <- data.table(
        chr=seqGetData(genofile, "chromosome"),
        pos=seqGetData(genofile, "position"),
        nAlleles=seqGetData(genofile, "$num_allele"),
        id=seqGetData(genofile, "variant.id")) %>%
    mutate(SNP_id = paste(chr, pos, sep = "_"))

#--------------------------------------------------------------------------------

# Extract annotation data for each SNP of interest (note: did not reannotate SNPs after added SNPs on either end of window)

annotation <- foreach(i=1:dim(outlier.win.SNPs)[1], .combine = "rbind", .errorhandling = "remove")%do%{

  message(i)
  # Reset filter
  seqResetFilter(genofile)
  # Extract SNP_id for SNP i
  tmp.i = outlier.win.SNPs[i,]$SNP_id
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

# Join annotation and SNP information
outlier.win.SNPs.annotated <- left_join(outlier.win.SNPs, annotation, by = join_by(SNP_id), relationship = "many-to-many")
outlier.win.SNPs.annotated <- outlier.win.SNPs.annotated %>% distinct()
# Write output
write.csv(outlier.win.SNPs.annotated, "data/processed/baypass/window_summary/outlier.win.SNPs.annotated.ph_mean.csv", row.names = F, quote = F)

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# Load windows
win.out.order.outliers <- read.csv("data/processed/baypass/window_summary/window_100kb_analysis_ph_mean_outliers.csv", header=T)

# Extract just top window (ntLink_3821, pos: 268 to 100162)
top.win <- win.out.order.outliers[which(win.out.order.outliers$rnp.binom.POD == min(win.out.order.outliers$rnp.binom.POD)),]

# Extract bf for top win
outlier.win.SNPs.top.win <- outlier.win.SNPs[which(outlier.win.SNPs$chr == top.win$chr & 
      outlier.win.SNPs$pos >= top.win$pos_min & 
      outlier.win.SNPs$pos <= top.win$pos_max), ]

# Extract SNPs within top outlier window (i.e., on chr ntLink_3821) - faster way
#outlier.win.ntlink3821 <- outlier.win.SNPs %>% filter(chr == "ntLink_3821")

# Graph BF for each SNP
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_topwin.pdf", width = 8.25, height = 4.5)
ggplot(outlier.win.SNPs.top.win, aes(y=bf_db.mean, x=pos/1000)) + labs(x="Position (kb)", y="BF") +
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.8, size=3.5, aes(colour = cut(bf_db.mean, c(-Inf, 20, Inf)))) + 
  scale_color_manual(values = c("(-Inf,20]" = "black", "(20, Inf]" = "blue")) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_topwin_wider.pdf", width = 10, height = 4.5)
ggplot(outlier.win.SNPs.top.win, aes(y=bf_db.mean, x=pos/1000)) + labs(x="Position (kb)", y="BF") +
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.8, size=3.5, aes(colour = cut(bf_db.mean, c(-Inf, 20, Inf)))) + 
  scale_color_manual(values = c("(-Inf,20]" = "black", "(20, Inf]" = "blue")) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()








