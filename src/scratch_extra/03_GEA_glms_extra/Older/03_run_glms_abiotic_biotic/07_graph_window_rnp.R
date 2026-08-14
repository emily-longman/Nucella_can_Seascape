# Graph window rnp

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
out_fig_dir <- paste("output/figures/GEA/glms/glms_window_summary")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load and merge data

# Create list of file names
path <- paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic_FET/")
file_names = as.list(dir(path = path, pattern = "glm_window_chunks_*"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0(paste("data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_abiotic_biotic_FET/"), x))))

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
save(win.out.order, file = "data/processed/GEA/glms/glms_window_summary/window_analysis_abiotic_biotic.RData")
load("data/processed/GEA/glms/glms_window_summary/window_analysis_abiotic_biotic.RData")

# ================================================================================== #

# Create unique Chromosome number
win.out.order.chr.unique <- unique(win.out.order$chr)
win.out.order$chr.unique <- as.numeric(factor(win.out.order$chr, levels = win.out.order.chr.unique))

# Graph rnp p
#pr.i = 0.01
pr.i = 0.05

# Graph rnp geompoint
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_abiotic_biotic_geompoint.pdf", width = 14, height = 6)
ggplot(win.out.order, aes(y=-log10(p.fet), x=chr.unique)) + 
  geom_point(alpha=0.8, size=1.6, aes(col = meanxtx)) +
  #scale_color_gradient(low = "grey", high = "blue") +
  scale_color_gradientn(colours=brewer.pal(9, "Greys")) +
  geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=26) #+ theme(legend.position = "none")
dev.off()

# Graph rnp geomline
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_abiotic_biotic_geomline.pdf", width = 12, height = 3)
ggplot(win.out.order, aes(y=-log10(p.fet), x=chr.unique)) + xlab("Scaffold") +
  geom_line( ) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=18) + theme(legend.position = "none") + theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()

# Graph Windows
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_abiotic_biotic_geomline_bywindows.pdf", width = 12, height = 3)
ggplot(win.out.order, aes(y=-log10(p.fet), x=win)) + xlab("Window (100kb)") +
  geom_line( ) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed") +
  theme_bw(base_size=18) + theme(legend.position = "none") #+ theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()
pdf("output/figures/GEA/glms/glms_window_summary/glm_window_abiotic_biotic_geomline_bywindows_alt.pdf", width = 20, height = 4)
ggplot(win.out.order, aes(y=-log10(p.fet), x=win)) + xlab("Window (100kb)") +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,2500, 5000, 7500, 10000, 12500)) + 
  #scale_y_continuous(limits=c(0,121), breaks=c(0, 30, 60, 90, 120)) +
  geom_line(linewidth=2) + geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed", linewidth=1.5) +
  theme_linedraw(base_size=28) + theme(legend.position = "none") #+ theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()

pdf("output/figures/GEA/glms/glms_window_summary/glm_window_abiotic_biotic_geompoint_bywindows.pdf", width = 20, height = 4)
ggplot(win.out.order, aes(y=-log10(p.fet), x=win)) + xlab("Window (100kb)") +
  geom_point(alpha=0.8, size=2, aes(col = meanxtx)) +
  #scale_color_gradient(low = "grey", high = "blue") +
  scale_color_gradientn(colours=brewer.pal(9, "Greys")) +
  scale_x_continuous(expand = c(0, 0), breaks=c(0,2500, 5000, 7500, 10000, 12500)) + 
  #scale_y_continuous(limits=c(0,121), breaks=c(0, 30, 60, 90, 120)) +
  geom_hline(yintercept=-log10(pr.i), col="red", linetype="dashed", linewidth=1.5) +
  theme_linedraw(base_size=28) #+ theme(legend.position = "none") #+ theme(plot.margin = margin(t = 40, r = 30, b = 20, l = 20, unit = "pt"))
dev.off()

# ================================================================================== #

# Extract outliers
#win.out.order.outliers <- win.out.order %>% filter(-log10(rnp.binom.p) > -log10(pr.i))
# Extract top outliers - those with enrichment >80
win.out.order.outliers <- win.out.order %>% filter(-log10(p.fet) > 1.5)

# ================================================================================== #

# Load data
load("data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/glm.output.all.Rdata")

# Identify SNPs in the outlier windows -> 45,095 SNPs
outlier.win.SNPs <- foreach(win.i=unique(win.out.order.outliers$win), .combine="rbind", .errorhandling="remove")%do%{
    
    # Extract window
    win.tmp <- win.out.order.outliers[which(win.out.order.outliers$win==win.i),]

    # Extract SNPs in window from baypass output
    glm.tmp <- o.all %>% filter(
        o.all$chr == win.tmp$chr & 
        o.all$pos >= win.tmp$pos_min &
        o.all$pos <= win.tmp$pos_max)
}

# Add SNP_id column
outlier.win.SNPs <- outlier.win.SNPs %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

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

# ================================================================================== #

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
write.csv(outlier.win.SNPs.annotated, "data/processed/GEA/glms/glms_window_summary/outlier_win_SNPs_abiotic_biotic_annotated_FET.csv", row.names = F, quote = F)
