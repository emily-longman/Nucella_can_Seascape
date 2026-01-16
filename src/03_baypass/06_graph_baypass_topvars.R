# Graph window analysis

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

# Load POD threshold
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")

# ================================================================================== #
 
# Read in Baypass pH mean files

# Load replicate BF outpu
bf.ph.mean <- foreach(i=1:4, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run", i, "_summary_betai_reg.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Change column names
setnames(bf.ph.mean, "BF(dB)", "bf_db")

# ================================================================================== #

# Average BF across replicate runs
bf.ph.mean.sum <- bf.ph.mean %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(bf_db.mean = mean(bf_db), bf_db.median = median(bf_db),
            bf_db.var=var(bf_db), eBPis.mean=mean(eBPis), eBPis.median=median(eBPis), eBPis.var=var(eBPis))

# ================================================================================== #

# Graph BF of individual variables
pdf("output/figures/baypass/baypass_BF_ph_mean_rep.means.pdf", width = 8, height = 8)
ggplot(bf.ph.mean.sum[which(bf.ph.mean.sum$bf_db > 0),], aes(y=bf_db.mean, x=chr)) + 
  labs(x = "Position", y = "BF (in dB)") +
  geom_point(alpha=0.8) + 
  geom_hline(yintercept=median(bf.sim.thr[thr==thr.use]$bf_db.mean)) +
  theme_classic(base_size = 20) + 
  theme(panel.spacing = unit(0.5, "lines"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12)) 
dev.off()


# ================================================================================== #

# Read in XtX files and graph - abiotic
foreach(w=abiotic_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/abiotic/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_summary_pi_xtx.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_xtx <- read.table(full_path, header=T)

    # Add variable column
    baypass_xtx <- baypass_xtx %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_xtx <- cbind(snp.meta, baypass_xtx)

    # Graph XtX for var
    pdf(paste0("output/figures/baypass/baypass_xtx_", w,".pdf"), width = 10, height = 5)
    plot(baypass_xtx$XtXst)
    dev.off()

    # Graph outliers for var
    pdf(paste0("output/figures/baypass/baypass_xtx_outliers_", w,".pdf"), width = 10, height = 5)
    par(mar=c(5,5,4,1)+.1) #Adjust margins
    plot(baypass_xtx$log10.1.pval., ylab=expression(-log[10](italic(p))), xlab="Position")
    abline(h=-log10(0.001), lty=2, col="red") #0.001 p-value threshold
    abline(h=-log10(0.05/dim(baypass_xtx)[1]), lty=1, col="red") # Bonferroni threshold
    dev.off()

}

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# Biotic data - M californianus thk

# Read all files and add variable name
baypass.bf <- foreach(w=Mcali_var, .combine="rbind", .errorhandling = "remove")%do%{  
    
    # State which variable loading
    message(w)

    # Path to file
    path <- paste("data/processed/baypass/biotic/Mcali/", w, sep = "")
    # File name
    file_name <- dir(path = path, pattern = "*_betai.out")
    # Full path
    full_path <- file.path(path, file_name)
    # Load file
    baypass_betai <- read.table(full_path, header=T)

    # Add variable column
    baypass_betai <- baypass_betai %>% mutate(variable = w)

    # Join baypass results with snp metadata
    baypass_betai_pos <- cbind(snp.meta, baypass_betai)
}

# Change column names
setnames(bf.POD, "BF(dB)", "bf_db")

# ================================================================================== #

# Explore patterns in bayes factors

# pH mean
bf.ph.mean <- baypass.bf[which(baypass.bf$variable == "ph_mean"),]
bf.ph.mean.extreme <- bf.ph.mean[which(bf.ph.mean$bf_db>50),]
bf.ph.mean.extreme <- bf.ph.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.ph.mean.extreme, "data/processed/baypass/bf.ph.mean.extremes.csv", row.names = F, quote = F)

# Temp mean
bf.thetao.mean <- baypass.bf[which(baypass.bf$variable == "thetao_mean"),]
bf.thetao.mean.extreme <- bf.thetao.mean[which(bf.thetao.mean$bf_db>50),]
bf.thetao.mean.extreme <- bf.thetao.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.thetao.mean.extreme, "data/processed/baypass/bf.thetao.mean.extremes.csv", row.names = F, quote = F)

# Temp min
bf.thetao.min <- baypass.bf[which(baypass.bf$variable == "thetao_min"),]
bf.thetao.min.extreme <- bf.thetao.min[which(bf.thetao.min$bf_db>50),]
bf.thetao.min.extreme <- bf.thetao.min.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.thetao.min.extreme, "data/processed/baypass/bf.thetao.min.extremes.csv", row.names = F, quote = F)

# Chl mean
bf.chl.mean <- baypass.bf[which(baypass.bf$variable == "chl_mean"),]
bf.chl.mean.extreme <- bf.chl.mean[which(bf.chl.mean$bf_db>50),]
bf.chl.mean.extreme <- bf.chl.mean.extreme %>% mutate(SNP_id = paste(chr, pos, sep = "_"))
write.csv(bf.chl.mean.extreme, "data/processed/baypass/bf.chl.mean.extreme.csv", row.names = F, quote = F)

# Note: All have a max Bayes Factor value of 52.96447989

# No overlap between the SNPs in the ph mean and thetao mean nor
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.thetao.mean.extreme$SNP_id)
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.chl.mean.extreme$SNP_id)
bf.ph.mean.extreme.overlap <- bf.ph.mean.extreme %>% filter(SNP_id %in% bf.thetao.min.extreme$SNP_id)


# ================================================================================== #
# ================================================================================== #

