# Graph Baypass XtX output.

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
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen', 'tidyverse', 'foreach', 'WriteXLS'))
library(data.table)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(geigen)
library(tidyverse)
library(foreach)
library(WriteXLS)

# ================================================================================== #

# Read in Baypass results
NC.omega <- as.matrix(read.table("data/processed/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out"))
XtX <- read.table("data/processed/outlier_analyses/baypass/xtx/NC_baypass_core_summary_pi_xtx.out", header=T)
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
snp.meta <- read.table("data/processed/outlier_analyses/baypass/snpdet", header=F)

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #

# Check the behavior of the p-values associated to the XtXst estimator
pdf("output/figures/outlier_analyses/Baypass_xtx_hist.pdf", width = 5, height = 5)
hist(10**(-1*XtX$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/outlier_analyses/Baypass_xtx.pdf", width = 5, height = 5)
plot(XtX$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/outlier_analyses/Baypass_xtx_outliers.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="-log10(XtX P-value)")
abline(h=3, lty=2, col="red") #0.001 p-value threshold
dev.off()

# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Identify bonferroni outliers

# Graph Bonferroni outliers (0.05 divided by the number of SNPs tested)
pdf("output/figures/outlier_analyses/Baypass_xtx_outliers_bonferroni.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="-log10(XtX P-value)" , xlab="Position")
abline(h=-log10(0.05/dim(XtX)[1]), col="red")
dev.off()

# Identify bonferroni significant SNPs -- 8 SNPS
bonf.sig.SNPs <- SNP.XtX.dt[which(SNP.XtX.dt$log10.1.pval. >= -log10(0.05/dim(SNP.XtX.dt)[1])),]

# Write file of bonferroni outlier SNPs
write.csv(bonf.sig.SNPs, "data/processed/outlier_analyses/baypass/Outlier_SNPs/Nucella_outlier_SNPs_bonferroni.csv", row.names=FALSE)

# ================================================================================== #

# Rank normalize p-values

L = dim(SNP.XtX.dt)[1]

# Calculate rank    
SNP.XtX.dt %>%
arrange(-log10.1.pval.) %>% 
as.data.frame() %>%
mutate(rank = seq(from=1,  to = L, by = 1)) %>% 
mutate(rank_norm = rank/L) %>% 
group_by(chr) %>%
arrange(as.numeric(pos))->
inner.rnf

# Window and step size
win.bp = 100000
step.bp = 50000

# Generate windows
wins <- foreach(chr.i=unique(SNP.XtX.dt$chr),
                .combine="rbind", 
                .errorhandling="remove")%do%{
                  
                  message(chr.i)
                  tmp <-  inner.rnf %>%
                    filter(chr == chr.i)

                    if(max(tmp$pos) <= step.bp){
                      data.table(chr=chr.i,
                             start=min(tmp$pos),
                             end=max(tmp$pos))  
                    } else if(max(tmp$pos) > step.bp){
                  
                  data.table(chr=chr.i,
                             start=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp),
                             end=seq(from=min(tmp$pos), to=max(tmp$pos)-win.bp, by=step.bp) + win.bp)
                } }

wins[,i:=1:dim(wins)[1]]

# Save windows
save(wins, file="data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_windows.RData")
# Reload windows
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_windows.RData")

# ================================================================================== #

# Start the summarization process

# Set p=0.01
win.out <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- inner.rnf %>%
    filter(chr == wins[win.i]$chr) %>%
    filter(pos >= wins[win.i]$start & pos <= wins[win.i]$end)
  
  pr.i <- c(0.01)
  
  win.tmp %>% 
    filter(!is.na(rank_norm)) %>%
    summarise(chr = wins[win.i]$chr,
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              win=win.i,
              pr=pr.i,
              rnp.pr=c(mean(rank_norm<=pr.i)),
              rnp.binom.p=c(binom.test(sum(rank_norm<=pr.i), 
                                       length(rank_norm), pr.i)$p.value),
              max.p=max(log10.1.pval.),
              nSNPs = n(),
              sum.rnp=sum(rank_norm<=pr.i),
    )  -> win.out
}

# Save window out
save(win.out, file="data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_win.out.RData")
# Reload windows
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_win.out.RData")

###################################################################################

# Set p=0.001
win.out.001 <- foreach(win.i=1:dim(wins)[1], 
                   .errorhandling = "remove",
                   .combine = "rbind"
)%do%{
  
  message(paste(win.i, dim(wins)[1], sep=" / "))
  
  
  win.tmp <- inner.rnf %>%
    filter(chr == wins[win.i]$chr) %>%
    filter(pos >= wins[win.i]$start & pos <= wins[win.i]$end)
  
  pr.i <- c(0.001)
  
  win.tmp %>% 
    filter(!is.na(rank_norm)) %>%
    summarise(chr = wins[win.i]$chr,
              pos_mean = mean(pos),
              pos_min = min(pos),
              pos_max = max(pos),
              win=win.i,
              pr=pr.i,
              rnp.pr=c(mean(rank_norm<=pr.i)),
              rnp.binom.p=c(binom.test(sum(rank_norm<=pr.i), 
                                       length(rank_norm), pr.i)$p.value),
              max.p=max(log10.1.pval.),
              nSNPs = n(),
              sum.rnp=sum(rank_norm<=pr.i),
    )  -> win.out.001
}

# Save window out
save(win.out.001, file="data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_win.out.001.RData")
# Reload windows
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/baypass_win.out.001.RData")

# ================================================================================== #

# pval=0.01

# Graph 

# Create unique Chromosome number 
chr.unique <- unique(win.out$chr)
win.out$chr.unique <- as.numeric(factor(win.out$chr, levels = chr.unique))

# Graph based on chr
pdf("output/figures/outlier_analyses/Baypass_rnp.pdf", width = 5, height = 5)
ggplot(win.out, aes(y=-log10(rnp.binom.p), x=chr.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()
# Graph based on chr - filter for windows with >50 SNPs and max p for that window > 2
pdf("output/figures/outlier_analyses/Baypass_rnp_filt.pdf", width = 5, height = 5)
ggplot(filter(win.out, nSNPs > 50 & max.p > 2), aes(y=-log10(rnp.binom.p), 
x=chr.unique)) + 
  geom_point(alpha=0.8) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()
# Graph based on chr - filter for windows with >50 SNPs and max p for that window > 2 (size of dot is max p for contig)
pdf("output/figures/outlier_analyses/Baypass_rnp_filt_size_max.p.pdf", width = 5, height = 5)
ggplot(filter(win.out, nSNPs > 50 & max.p > 2), aes(y=-log10(rnp.binom.p), 
x=chr.unique, size = max.p)) + 
  geom_point(alpha=0.8) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()
# Graph based on position
pdf("output/figures/outlier_analyses/Baypass_rnp_pos.pdf", width = 5, height = 5)
ggplot(filter(win.out, nSNPs > 50 & max.p > 2), 
aes(y=-log10(rnp.binom.p), x=pos_mean/1e6, size = max.p)) + 
  geom_point( alpha=0.8, 
  #size=1.3
  ) + 
  geom_hline(yintercept = -log10(0.01), color="red") +
  theme_bw()
dev.off()

# ================================================================================== #

# Identify contigs with highly significant rnp p
win.out.sig <- win.out[which(-log10(win.out$rnp.binom.p) >= -log10(0.01)),]

# Create outlier SNP list (total of 838,918 SNPs)
SNPs.Interest <- foreach(i=1:dim(win.out.sig)[1], .combine = "rbind")%do%{
  tmp.snps <- inner.rnf %>%
  filter(chr == win.out.sig[i,]$chr) %>%
  filter(pos >= win.out.sig[i,]$pos_min & pos <= win.out.sig[i,]$pos_max)
}

# Write file of outlier SNPs
write.csv(SNPs.Interest, "data/processed/outlier_analyses/baypass/Outlier_SNPs/Nucella_outlier_SNPs.csv", row.names=FALSE)

# Save SNPs of interest
save(SNPs.Interest, file="data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.01.RData")
# Reload SNPs of interest
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.01.RData")

# ================================================================================== #
# ================================================================================== #
# ================================================================================== #

# pval=0.001

# Graph 

# Create unique Chromosome number 
chr.unique <- unique(win.out.001$chr)
win.out.001$chr.unique <- as.numeric(factor(win.out.001$chr, levels = chr.unique))

# Graph based on chr
pdf("output/figures/outlier_analyses/Baypass_rnp_pval_0.001.pdf", width = 5, height = 5)
ggplot(win.out.001, aes(y=-log10(rnp.binom.p), x=chr.unique)) + 
  geom_point(col="black", alpha=0.8, size=1.3) + 
  geom_hline(yintercept = -log10(0.001), color="red") +
  theme_bw()
dev.off()
# Graph based on chr - filter for windows with >50 SNPs and max p for that window > 2
pdf("output/figures/outlier_analyses/Baypass_rnp_filt_pval_0.001.pdf", width = 5, height = 5)
ggplot(filter(win.out.001, nSNPs > 50 & max.p > 2), aes(y=-log10(rnp.binom.p), 
x=chr.unique, color = max.p)) + 
  geom_point(alpha=0.8) + 
  geom_hline(yintercept = -log10(0.001), color="red") +
  theme_bw()
dev.off()

# ================================================================================== #

# Identify contigs with highly significant rnp p

# Identity outlier SNPs 178 SNPS
win.out.001.sig <- win.out.001[which(-log10(win.out.001$rnp.binom.p) >= -log10(0.001)),]

# Create contig list SNP list  (total of 70,206 SNPs)
SNPs.Interest.pval.001 <- foreach(i=1:dim(win.out.001.sig)[1], .combine = "rbind")%do%{
  tmp.snps <- inner.rnf %>%
  filter(chr == win.out.001.sig[i,]$chr) %>%
  filter(pos >= win.out.001.sig[i,]$pos_min & pos <= win.out.001.sig[i,]$pos_max)
}

# Write file of outlier SNPs
write.csv(SNPs.Interest.pval.001, "data/processed/outlier_analyses/baypass/Outlier_SNPs/Nucella_outlier_SNPs_pval.0.001.csv", row.names=FALSE)

# Save SNPs of interest
save(SNPs.Interest.pval.001, file="data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.001.RData")
# Reload SNPs of interest
load("data/processed/outlier_analyses/baypass/Outlier_SNPs/SNPs.Interest.pval.001.RData")
