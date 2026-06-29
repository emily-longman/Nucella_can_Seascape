# Graph pH outliers

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'poolfstat', 'RColorBrewer', 'viridis'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(poolfstat)
library(RColorBrewer)
library(viridis)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass/outliers")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load mean bf data from 5 baypass runs
load("data/processed/baypass/abiotic/bf.ph.mean.sum.Rdata")
# Load windows
win.out.order.outliers <- read.csv("data/processed/baypass/window_summary/window_100kb_analysis_ph_mean_outliers.csv", header=T)
# Extract just top window (ntLink_3821, pos: 268 to 100162)
top.win <- win.out.order.outliers[which(win.out.order.outliers$rnp.binom.POD == min(win.out.order.outliers$rnp.binom.POD)),]

# Load POD thresholds
load("data/processed/baypass/abiotic/ph_mean_POD_thr.Rdata")

# Load SNPs of interest
#bf.ph.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header=T)
# Extract just those on g27343
#bf.ph.mean.g27343 <- bf.ph.mean.sum.outliers.annotated %>% filter(Gene_Name == "g27343")

# Load bio-oracle environmental data
bio_oracle_sites_2010 <- read.csv("data/processed/GEA/enviro_data/Bio-oracle/bio_oracle_sites_2010.csv", header=T)
# Extract just mean pH and rename location to Site
ph <- bio_oracle_sites_2010[,c(1,2,3,11,13)]
ph <- ph %>% rename(Site = location)

# ================================================================================== #

# Extract just BF for top win
bf.ph.mean.sum.top.win <- bf.ph.mean.sum[which(bf.ph.mean.sum$chr == top.win$chr & 
      bf.ph.mean.sum$pos >= top.win$pos_min & 
      bf.ph.mean.sum$pos <= top.win$pos_max), ]

bf.ph.mean.sum.top.win <- bf.ph.mean.sum.top.win %>% mutate(SNP_id = paste(chr, pos, sep = "_"))

# ================================================================================== #

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Subset pooldata for Baypass outlier SNPs
# Extract SNP info for all SNPs and make snp_id column
pooldata@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info
# Rename columns
names(snp.info)[1:2] = c("chr","pos")
# Make snp_id column
snp.info %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info

# Filter pooldata for Baypass SNPs of interest
selected_SNPs_ph <- snp.info %>% filter(SNP_id %in% bf.ph.mean.sum.top.win$SNP_id)
# Get index of SNPs
selected_SNPs_ph_index <- as.integer(selected_SNPs_ph$rs.id)
# Subset the pooldata object using the selected SNP indices
pooldata_ph <- pooldata.subset(pooldata, snp.index = selected_SNPs_ph_index)

# Extract and manipulate snp info for significant SNPs

# Extract SNP info for significant SNPs
pooldata_ph@snp.info %>%
  as.data.frame() %>% mutate(rs.id = rownames(.)) -> snp.info.ph
# Rename columns
names(snp.info.ph)[1:2] = c("chr","pos")
# Make snp_id column
snp.info.ph %>% mutate(SNP_id = paste(chr, pos, sep = "_")) -> snp.info.ph

# ================================================================================== #

# Extract and manipulate coutn and coverage for significant SNPs

# Extract read count data for SNPs
ref_count <- pooldata_ph@refallele.readcount
ref_count %>% as.data.frame -> count
names(count) = c(pooldata_ph@poolnames)
count$SNP_id <- snp.info.ph$SNP_id
count.melt <- reshape2::melt(count, id = "SNP_id", variable.name = "Site", value.name = "Count")

# Extract coverage data for SNPs
coverage <- pooldata_ph@readcoverage
coverage %>% as.data.frame -> cov
names(cov) = c(pooldata_ph@poolnames)
cov$SNP_id <- snp.info.ph$SNP_id
cov.melt <- reshape2::melt(cov, id = "SNP_id", variable.name = "Site", value.name = "Cov")

# Calculate allele frequency for SNPs
allele_freqs <- ref_count/coverage
# Change to data frame
allele_freqs %>% as.data.frame -> afs
# Rename columns (19 sites)
names(afs) = c(pooldata_ph@poolnames)
# Add SNP_id
afs$SNP_id <- snp.info.ph$SNP_id
# Change format
afs.melt <- reshape2::melt(afs, id = "SNP_id", variable.name = "Site", value.name = "AF")

# Join
ph.all <- left_join(count.melt, cov.melt)
ph.all <- left_join(ph.all, afs.melt)

# Look at table for one site (top missense mutation)
ph.all[which(ph.all$SNP_id=="ntLink_3821_41192"),]

# Join with pH data
afs.ph <- left_join(ph.all, ph, by="Site")


# Make Site an ordered factor
lat.order <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
afs.ph <- afs.ph %>% mutate(Site = factor(Site, levels = lat.order))
# Add column that highlights N, and S
afs.ph <- afs.ph %>% mutate(shape = case_when(Site %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   Site %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))

# Color palette
viridiscolors <- viridis(n=19)

# Subset so only SNPs with BF > 20 or > 18, etc
afs.ph.BF22 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.sum.top.win[which(bf.ph.mean.sum.top.win$bf_db.mean > 22),]$SNP_id)
afs.ph.BF20 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.sum.top.win[which(bf.ph.mean.sum.top.win$bf_db.mean > 20),]$SNP_id)
afs.ph.BF19 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.sum.top.win[which(bf.ph.mean.sum.top.win$bf_db.mean > 19),]$SNP_id)
afs.ph.BF18 <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.sum.top.win[which(bf.ph.mean.sum.top.win$bf_db.mean > 18),]$SNP_id)

afs.ph.g27343.BF.POD <- afs.ph %>% filter(SNP_id %in% bf.ph.mean.sum.top.win[which(bf.ph.mean.sum.top.win$bf_db.mean > bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)]),]$SNP_id)
# Save
save(afs.ph.g27343.BF.POD, file = "data/processed/baypass/afs.ph.g27343.BF.POD.RData")


# Graph and color by Site - BF 22
pdf("output/figures/baypass/outliers/pH_g27343_BF22.pdf", width = 9, height = 3.75)
ggplot(afs.ph.BF22, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.8, size = 6) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Frequency", y="mean pH") +
  facet_wrap(~SNP_id, ncol = 4) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) + scale_y_continuous(breaks = c(7.95, 8.00)) +
  theme_linedraw(base_size = 30) + theme(legend.position="none") + 
  theme(strip.background = element_rect(fill = "White"), strip.text = element_text(color = "black"))
dev.off()
# Graph and color by Site - BF20
pdf("output/figures/baypass/outliers/pH_g27343_BF20.pdf", width = 18, height = 4)
ggplot(afs.ph.BF20, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.7, size = 5) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Frequency", y="mean pH") +
  facet_wrap(~SNP_id, ncol = 4) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) +
  theme_bw(base_size = 26) + theme(legend.position="none")
dev.off()
pdf("output/figures/baypass/outliers/pH_g27343_BF20_wider.pdf", width = 9.435, height = 6.5)
ggplot(afs.ph.BF20, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.7, size = 7) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Frequency", y="mean pH") +
  facet_wrap(~SNP_id, ncol = 2) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) + scale_y_continuous(limits = c(7.92, 8.03), breaks = c(7.95, 8.0)) +
  theme_linedraw(base_size = 30) + theme(strip.background =element_rect(fill="grey"), strip.text = element_text(colour = 'black')) +
  theme(legend.position="none")
dev.off()
pdf("output/figures/baypass/outliers/pH_g27343_BF20_taller.pdf", width = 4.5, height = 8.5)
ggplot(afs.ph.BF20, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.7, size = 7) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Frequency", y="mean pH") +
  facet_wrap(~SNP_id, ncol = 1) + scale_fill_manual(values = viridiscolors) +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) + scale_y_continuous(limits = c(7.92, 8.03), breaks = c(7.95, 8.0)) +
  theme_linedraw(base_size = 30) + 
  theme(strip.background = element_rect(fill="grey"), strip.text = element_text(colour = 'black', size = 18, margin = margin(t=6, r=6, b=6, l=6))) +
  theme(legend.position="none")
dev.off()
# Graph and color by Site - BF19
pdf("output/figures/baypass/outliers/pH_g27343_BF19.pdf", width = 16, height = 9)
ggplot(afs.ph.BF19, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.7, size = 6) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Freq", y="mean pH") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) + scale_y_continuous(breaks = c(7.95, 8.00)) +
  facet_wrap(~SNP_id) + scale_fill_manual(values = viridiscolors) +
  theme_bw(base_size = 26)
dev.off()
# Graph and color by Site - BF18
pdf("output/figures/baypass/outliers/pH_g27343_BF18.pdf", width = 18, height = 12)
ggplot(afs.ph.BF18, aes(x=AF, y=ph_mean, shape=shape, fill=Site)) +
  geom_point(alpha=0.7, size = 6) + scale_shape_manual(values = c(21, 23)) + labs(x="Allele Freq", y="mean pH") +
  scale_x_continuous(breaks = seq(0, 1, by = 0.5)) + scale_y_continuous(breaks = c(7.95, 8.00)) +
  facet_wrap(~SNP_id) + scale_fill_manual(values = viridiscolors) +
  theme_bw(base_size = 26) + guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Graph BF for top win

# Graph BF for each SNP within window
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_topwin.pdf", width = 8, height = 4)
ggplot(bf.ph.mean.sum.top.win, aes(y=bf_db.mean, x=pos/1000)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.75, size=3.5) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_topwin_wider.pdf", width = 10, height = 4)
ggplot(bf.ph.mean.sum.top.win, aes(y=bf_db.mean, x=pos/1000)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.75, size=3.5) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_topwin_wider_coltopSNPs.pdf", width = 10, height = 4)
ggplot(bf.ph.mean.sum.top.win, aes(y=bf_db.mean, x=pos/1000)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.75, size=4.5, shape=21, aes(fill = cut(bf_db.mean, c(-Inf, 22, Inf)))) + 
  scale_fill_manual(values = c("(-Inf,22]" = "black", "(22, Inf]" = "#24dafa")) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Load and average xtx

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load xtx output for 5 replicate Baypass runs
baypass.ph.xtx <- foreach(i=1:5, .combine = rbind)%do%{
    message(i)
    tmp <- fread(paste("data/processed/baypass/abiotic/ph_mean/NC_abiotic_ph_mean_run", i, "_summary_pi_xtx.out", sep=""))
    tmp[,rep:=i]
    tmp <- cbind(snp.meta, tmp)
    return(tmp)
}

# Rename p val
baypass.ph.xtx <- baypass.ph.xtx %>% rename(log10.1.pval. = "log10(1/pval)")

# Average across replicate runs
baypass.ph.xtx.sum <- baypass.ph.xtx %>% group_by(chr, pos, allele1, allele2, MRK) %>% 
    reframe(M_P_mean = mean(M_P), SD_P_mean = mean(SD_P), M_XtX_mean = mean(M_XtX), 
    SD_XtX_mean = mean(SD_XtX), XtXst_mean = mean(XtXst), log10.1.pval_mean = mean(log10.1.pval.))

# Save
save(baypass.ph.xtx.sum, file = "data/processed/baypass/abiotic/baypass.ph.xtx.sum.RData")
load("data/processed/baypass/abiotic/baypass.ph.xtx.sum.RData")

####

# Load PODs and average

# Create list of file names
file_names = as.list(dir(path = 'data/processed/baypass/abiotic/ph_mean_POD/', pattern = "*summary_pi_xtx.out"))
file_names_v = as.vector(unlist(lapply(file_names, function(x) paste0('data/processed/baypass/abiotic/ph_mean_POD/', x))))

# Read all the files and add a column with the run
baypass.ph.xtx.POD <- foreach(w=file_names_v, .combine = rbind)%do%{  
    # State which file loading
    message(w)
    # Load file
    tmp = fread(w, header=T)
    # Add column with identifier
    tmp <- tmp %>% mutate(run = w) %>% mutate(run = str_remove(run, pattern = "data/processed/baypass/abiotic/ph_mean_POD/NC_abiotic_ph_mean_POD_run*"))
    # Remove end of chunk name
    tmp <- tmp %>% mutate(run = str_remove(run, pattern = "_summary_pi_xtx.out"))
    #Return
    return(tmp)
}

# Calculate quantiles for each POD
baypass.ph.xtx.POD.sum <- baypass.ph.xtx.POD %>% group_by(run) %>% reframe(XtXst = quantile(XtXst, c(.95, .99, .999)), M_XtX = quantile(M_XtX, c(.95, .99, .999)), thr = c(.95, .99, .999)) %>% as.data.frame()

# Average quantiles across POD runs
baypass.ph.xtx.POD.thr <- baypass.ph.xtx.POD.sum %>% group_by(thr) %>% summarize(XtXst_mean=mean(XtXst), M_XtX_mean = mean(M_XtX))

# ================================================================================== #

# Graph xtx for top window

# Extract xtx for just top window
baypass.ph.xtx.sum.top.win <- baypass.ph.xtx.sum[which(baypass.ph.xtx.sum$chr == top.win$chr & 
      baypass.ph.xtx.sum$pos >= top.win$pos_min & 
      baypass.ph.xtx.sum$pos <= top.win$pos_max), ]

# Graph xtx for each SNP within window
pdf("output/figures/baypass/window_summary/baypass_ph_mean_xtx_topwin.pdf", width = 8, height = 4)
ggplot(baypass.ph.xtx.sum.top.win, aes(y=XtXst_mean, x=pos/1000)) + labs(x="Position (kb)", y=expression(italic("XtX")))+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.75, size=3.5) + ylim(0,50) +
  geom_hline(yintercept=baypass.ph.xtx.POD.thr$XtXst_mean[which(baypass.ph.xtx.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()
pdf("output/figures/baypass/window_summary/baypass_ph_mean_xtx_topwin_wider.pdf", width = 9.80595, height = 4)
ggplot(baypass.ph.xtx.sum.top.win, aes(y=XtXst_mean, x=pos/1000)) + labs(x="Position (kb)", y="XtX")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.75, size=4, shape=21, fill="black") + scale_y_continuous(limits=c(2,48), breaks=c(0,15,30,45)) +
  geom_hline(yintercept=baypass.ph.xtx.POD.thr$XtXst_mean[which(baypass.ph.xtx.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()

# Graph xtx for each SNP within window - geom_line
pdf("output/figures/baypass/window_summary/baypass_ph_mean_xtx_topwin_geomline.pdf", width = 8, height = 4.5)
ggplot(baypass.ph.xtx.sum.top.win, aes(y=XtXst_mean, x=pos/1000)) + labs(x="Position (kb)", y="XtX") +
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_line() + 
  geom_hline(yintercept=baypass.ph.xtx.POD.thr$XtXst_mean[which(baypass.ph.xtx.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()

# Graph CORRECTED xtx for each SNP within window
pdf("output/figures/baypass/window_summary/baypass_ph_mean_xtx_corrected_topwin.pdf", width = 8, height = 4.5)
ggplot(baypass.ph.xtx.sum.top.win, aes(y=M_XtX_mean, x=pos/1000)) + labs(x="Position (kb)", y="XtX corrected")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), fill="grey", alpha=0.5) +
  geom_point(alpha=0.8, size=3.5) + ylim(0,35) +
  geom_hline(yintercept=baypass.ph.xtx.POD.thr$M_XtX_mean[which(baypass.ph.xtx.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  theme_bw(base_size=30) + theme(legend.position = "none")
dev.off()





# Join BF and xtx
baypass.ph.top.win <- left_join(bf.ph.mean.sum.top.win, baypass.ph.xtx.sum.top.win)
# Graph
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_xtx_topwin_wider.pdf", width = 9.36, height = 4)
ggplot(baypass.ph.top.win, aes(y=bf_db.mean, x=pos/1000, col = M_XtX_mean)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), color="grey" , fill="grey", alpha=0.5) +
  #geom_point(alpha=0.75, size=6) + 
  geom_point(alpha=0.75, size=7, aes(shape = cut(bf_db.mean, c(-Inf, 20, Inf)))) + 
  #scale_color_gradient(low = "#d7dbf6", high = "blue", name=expression(italic("X"^T*"X"))) + 
  scale_color_gradientn(colours=rev(brewer.pal(9, "RdYlBu")), name=expression(italic("X"^T*"X"))) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()
pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_xtx_topwin_alt.pdf", width = 8.85, height = 4.5)
ggplot(baypass.ph.top.win, aes(y=bf_db.mean, x=pos/1000, col = M_XtX_mean)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), color="grey" , fill="grey", alpha=0.5) +
  #geom_point(alpha=0.75, size=6) + 
  geom_point(alpha=0.75, size=7, aes(shape = cut(bf_db.mean, c(-Inf, 20, Inf)))) + 
  #scale_color_gradient(low = "#d7dbf6", high = "blue", name=expression(italic("X"^T*"X"))) + 
  scale_color_gradientn(colours=rev(brewer.pal(9, "RdYlBu")), name=expression(italic("X"^T*"X"))) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()

pdf("output/figures/baypass/window_summary/baypass_ph_mean_BF_xtx_topwin_alt2.pdf", width = 8.85, height = 4.5)
ggplot(baypass.ph.top.win, aes(y=bf_db.mean, x=pos/1000, col = M_XtX_mean)) + labs(x="Position (kb)", y="BF")+
  geom_rect(aes(xmin=1/1000, xmax=58041/1000, ymin=-Inf, ymax=Inf), color="grey" , fill="grey", alpha=0.5) +
  #geom_point(alpha=0.75, size=6) + 
  geom_point(alpha=0.75, size=7, aes(shape = cut(bf_db.mean, c(-Inf, 20, Inf)))) + 
  geom_point(data = baypass.ph.top.win[which(baypass.ph.top.win$bf_db.mean>20),], size=7, aes(shape = cut(bf_db.mean, c(-Inf, 20, Inf)))) +
  geom_point(data = baypass.ph.top.win[which(baypass.ph.top.win$bf_db.mean>20),], shape = 24, size = 7, color = "black", fill = "transparent") +
  #geom_point(data = baypass.ph.top.win[which(baypass.ph.top.win$bf_db.mean>20),], aes(y=bf_db.mean, x=pos/1000, fill = M_XtX_mean), shape = 24, color = "black") +
  ylim(-13,26.5) +
  #scale_shape_manual(values = c(21, 24)) +
  #scale_color_gradient(low = "#d7dbf6", high = "blue", name=expression(italic("X"^T*"X"))) + 
  scale_color_gradientn(colours=rev(brewer.pal(9, "RdYlBu")), name=expression(italic("X"^T*"X"))) +
  scale_fill_gradientn(colours=rev(brewer.pal(9, "RdYlBu")), name=expression(italic("X"^T*"X"))) +
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], col="red", linetype="dashed") +
  geom_hline(yintercept=0, col="black", linetype="solid") +
  theme_linedraw(base_size=30) + theme(legend.position = "none")
dev.off()