# Graph xtx

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'ggplot2', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Load SeqArray
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install(version = "3.20")
#BiocManager::install("SeqArray")
library(SeqArray)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load Data

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# ================================================================================== #

# Load xtx output for 5 replicate Baypass runs
baypass.Mcali.xtx <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/biotic/Mcali_IntegratedThk/NC_biotic_Mcali_IntegratedThk_run", i, "_summary_pi_xtx.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Rename p val
baypass.Mcali.xtx <- baypass.Mcali.xtx %>% rename(log10.1.pval. = "log10(1/pval)")

# Average across replicate runs
baypass.Mcali.xtx.sum <- baypass.Mcali.xtx %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(M_P_mean = mean(M_P), SD_P_mean = mean(SD_P), M_XtX_mean = mean(M_XtX), 
    SD_XtX_mean = mean(SD_XtX), XtXst_mean = mean(XtXst), log10.1.pval_mean = mean(log10.1.pval.))

# Save
save(baypass.Mcali.xtx.sum, file = "data/processed/baypass/biotic/baypass.Mcali.xtx.sum.RData")
load("data/processed/baypass/biotic/baypass.Mcali.xtx.sum.RData")

####

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/biotic/Mcali_IntegratedThk_POD/', pattern = "*summary_pi_xtx.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/biotic/Mcali_IntegratedThk_POD/', x))))

# Read all the files and add a column with the run
baypass.Mcali.xtx.POD <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/baypass/biotic/Mcali_IntegratedThk_POD/NC_biotic_Mcali_IntegratedThk_POD_run*"))
    # Remove end of chunk name
    tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_pi_xtx.out"))
    #Return
    return(tmp)
}

# Calculate quantiles for each POD
baypass.Mcali.xtx.POD.sum <- baypass.Mcali.xtx.POD %>% group_by(run) %>% reframe(XtXst = quantile(XtXst, c(.95, .99, .999)), M_XtX = quantile(M_XtX, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()

# Average quantiles across POD runs
baypass.Mcali.xtx.POD.thr <- baypass.Mcali.xtx.POD.sum %>% group_by(thr) %>% summarize(XtXst_mean=mean(XtXst), M_XtX_mean = mean(M_XtX))


# Graph corrected xtx for each SNP within window
pdf("output/figures/baypass/baypass_Mcali_xtx.pdf", width = 12, height = 3)
ggplot(baypass.Mcali.xtx.sum, aes(y=M_XtX_mean, x=MRK/1000)) + labs(x="Position (kb)", y=expression(paste(italic("XtX"), "corrected")))+
  geom_point(alpha=0.8, size=3.5) + ylim(0,40) +
  geom_hline(yintercept=baypass.Mcali.xtx.POD.thr$M_XtX_mean[which(baypass.Mcali.xtx.POD.thr$thr==0.999)], col="#d900ff", linetype="dashed") +
  theme_bw(base_size=20) + theme(legend.position = "none")
dev.off()

# Extract outliers
baypass.Mcali.xtx.sum.outliers <- baypass.Mcali.xtx.sum[which(baypass.Mcali.xtx.sum$M_XtX_mean>30),]

# Create SNP_id
baypass.Mcali.xtx.sum.outliers <- baypass.Mcali.xtx.sum.outliers %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

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

annotation <- foreach(i=1:dim(baypass.Mcali.xtx.sum.outliers)[1], .combine = "rbind", .errorhandling = "remove")%do%{

  message(i)
  # Reset filter
  seqResetFilter(genofile)
  # Extract SNP_id for SNP i
  tmp.i = baypass.Mcali.xtx.sum.outliers[i,]$SNP_id
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
baypass.Mcali.xtx.sum.outliers.annotated <- left_join(baypass.Mcali.xtx.sum.outliers, annotation, by = join_by(SNP_id), relationship = "many-to-many")
baypass.Mcali.xtx.sum.outliers.annotated <- baypass.Mcali.xtx.sum.outliers.annotated %>% distinct()
# Write output
write.csv(baypass.Mcali.xtx.sum.outliers.annotated, "data/processed/baypass/baypass.Mcali.xtx.sum.outliers.annotated.csv", row.names = F, quote = F)
