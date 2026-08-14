# Calculate Fst statistics and gnerate PCA of all SNPs

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root_file(criterion = has_file("README.md")))
root_path <- find_root_file(criterion = has_file("README.md"))
# Set working directory as path from root
setwd(root_path)

# ================================================================================== #

# Load packages
install.packages(c('poolfstat', 'ggplot2', 'RColorBrewer'))
library(poolfstat)
library(ggplot2)
library(RColorBrewer)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 2, max.cov.per.pool = 200, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 16,398,082 SNPs for 19 Pools

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

#### 

# Try more stringent filters
#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf.gz", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 30, min.rc = 2, max.cov.per.pool = 100, min.maf = 0.1, nlines.per.readblock = 1e+06)
# Data consists of 2,671,913 SNPs for 19 Pools

#######
# Filtered 
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 2, max.cov.per.pool = 200, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 10,109,226 SNPs for 19 Pools

###
# LD pruned - can't seem to load (ERROR: No field containing allele depth (AD field) was detected in the vcf file)
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_LD/N.canaliculata_pops.plink.LDfiltered_0.8.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 30, min.rc = 2, max.cov.per.pool = 100, min.maf = 0.1, nlines.per.readblock = 1e+06)

# ================================================================================== #

# variant calling with bcftools 

# Read in data and filter
#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 25, min.rc = 2, max.cov.per.pool = 200, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 5,950,444 SNPs for 19 Pools

# ================================================================================== #

# Estimate genome wide Fst 

# Use computeFST function
pooldata.fst <- computeFST(pooldata,verbose=FALSE)
pooldata.fst$Fst 
# 0.5330188
# Relaxed: 0.5202263
# Stringent: 0.6536584
# bcftools: 0.5337656

# Block-Jackknife estimation of Fst standard error and confidence intervals
pooldata.fst.bjack <- computeFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)
pooldata.fst.bjack$Fst
#   Estimate  bjack mean  bjack s.e.     CI95inf     CI95sup 
# 0.533018756 0.533764540 0.001051527 0.531703547 0.535825534 
# Relaxed: 0.5202263006 0.5206428107 0.0009732689 0.5187352036 0.5225504178
# Stringent: 0.65365838 0.65880523 0.00487404 0.64925211 0.66835834 
# bcftools: 0.533765580 0.538158861 0.001886336 0.534461643 0.541856080 

# Compute multi-locus Fst over sliding window of SNPs
pooldata.fst.sliding.window <- computeFST(pooldata, sliding.window.size=100)

# Plot sliding window
pdf("output/figures/pop_structure/fst.sliding.window.pdf", width = 10, height = 10)
plot(pooldata.fst.sliding.window$sliding.windows.fvalues$CumMidPos/1e6, 
pooldata.fst.sliding.window$sliding.windows.fvalues$MultiLocusFst,
xlab="Cumulated Position (in Mb)", ylab="Multi-locus Fst", pch=16)
abline(h=pooldata.fst.sliding.window$Fst,lty=2, col="red") # Dashed red line indicates the estimated overall genome-wide Fst
dev.off()

# ================================================================================== #

# Estimate pairwise-population Fst

# Use compute.pairwiseFST function
pooldata.pairwisefst <- compute.pairwiseFST(pooldata, verbose=FALSE)

# Graph heatmap
pdf("output/figures/pop_structure/heatmap.pdf", width = 10, height = 10)
heatmap(pooldata.pairwisefst)
dev.off()

# Block-Jackknife estimation of pairwise Fst standard error and confidence intervals
pooldata.pairwisefst.bjack <- compute.pairwiseFST(pooldata, nsnp.per.bjack.block = 1000, verbose=FALSE)

# Estimated pairwise Fst are stored in the slot values: 5 first estimated pairwise
head(pooldata.pairwisefst.bjack@values)

# Graph estimated pairwise-population FST with their 95% confidence intervals 
pdf("output/figures/pop_structure/pairwise_Fst.pdf", width = 10, height = 10)
plot(pooldata.pairwisefst.bjack, cex=0.5)
dev.off()

# ================================================================================== #

# Principle Components Analysis with randomallele.pca

#PCA on the read count data (the object)
pooldata.pca = randomallele.pca(pooldata, main="Read Count data")

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.reorder <- mycolors[c(19,2,3,4,11,5,1,10,17,14,6,8,7,13,9,18,16,12,15)]

#colors.reorder.alphabetical <- mycolors[c(4, 11, 5, 1, 10, 17, 8, 18, 16, 12, 13, 6, 15, 14, 3, 2, 7, 19, 9)]

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/PCA_all_SNPs_PC1_PC2.pdf", width = 10, height = 10)
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 3, main="Read Count data")
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/PCA_all_SNPs_PC3_PC4.pdf", width = 10, height = 10)
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, main="Read Count data")
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plot names on PCAs

# Plotting PC1 and PC2
pdf("output/figures/pop_structure/PCA_all_SNPs_PC1_PC2.pdf", width = 10, height = 10)
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 3, main="Read Count data")
text(pooldata.pca$pop.loadings[,1], pooldata.pca$pop.loadings[,2], pooldata@poolnames)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

# Plotting PC3 and PC4
pdf("output/figures/pop_structure/PCA_all_SNPs_PC3_PC4.pdf", width = 10, height = 10)
pca <- plot(pooldata.pca$pop.loadings[,3],pooldata.pca$pop.loadings[,4],
xlab=paste0("PC",3," (",round(pooldata.pca$perc.var[3],2),"%)"),
ylab=paste0("PC",4," (",round(pooldata.pca$perc.var[4],2),"%)"),
col="black", bg=colors.reorder,pch=21, cex = 3, main="Read Count data")
text(pooldata.pca$pop.loadings[,3], pooldata.pca$pop.loadings[,4], pooldata@poolnames)
abline(h=0,lty=2,col="grey") ; abline(v=0,lty=2,col="grey")
dev.off()

## ggplot
pdf("output/figures/pop_structure/PCA_all_SNPs_PC1_PC2_ggplot.pdf", width = 10, height = 10)
ggplot(data=pooldata.pca, aes(x=pooldata.pca$pop.loadings[,1], y=pooldata.pca$pop.loadings[,2])) + 
geom_point(size=2, shape=21, color=colors.reorder) + 
xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) +
theme_classic()
dev.off()

