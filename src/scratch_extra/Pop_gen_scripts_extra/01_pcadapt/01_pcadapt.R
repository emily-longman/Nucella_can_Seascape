# Build pcadapt matrix 
# i.e., build a matrix of relative frequencies with n rows and L columns, 
# where n is the number of populations and L is the number of genetic markers

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
install.packages(c('poolfstat', 'pcadapt', 'RColorBrewer'))
library(poolfstat)
library(pcadapt)
library(RColorBrewer)
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
library(qvalue)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 26,206,958)
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 11,656,080 SNPs for 19 Pools

# ================================================================================== #

# S4 class to represent a pool-seq data set 
# Need to convert this to a matrix (follow: https://github.com/chaberko-lbbe/clec-poolseq)

# Specify reference allele for each population
ref_STR <- pooldata@refallele.readcount[,1]
ref_SLR <- pooldata@refallele.readcount[,2]
ref_SH <- pooldata@refallele.readcount[,3]
ref_ARA <- pooldata@refallele.readcount[,4]
ref_BMR <- pooldata@refallele.readcount[,5]
ref_CBL <- pooldata@refallele.readcount[,6]
ref_FC <- pooldata@refallele.readcount[,7]
ref_FR <- pooldata@refallele.readcount[,8]
ref_HZD <- pooldata@refallele.readcount[,9]
ref_SBR <- pooldata@refallele.readcount[,10]
ref_PSG <- pooldata@refallele.readcount[,11]
ref_KH <- pooldata@refallele.readcount[,12]
ref_STC <- pooldata@refallele.readcount[,13]
ref_PL <- pooldata@refallele.readcount[,14]
ref_VD <- pooldata@refallele.readcount[,15]
ref_OCT <- pooldata@refallele.readcount[,16]
ref_PB <- pooldata@refallele.readcount[,17]
ref_PGP <- pooldata@refallele.readcount[,18]
ref_PSN <- pooldata@refallele.readcount[,19]

# Calculate alternate allele for each population
alt_STR <- pooldata@readcoverage[,1] - pooldata@refallele.readcount[,1]
alt_SLR <- pooldata@readcoverage[,2] - pooldata@refallele.readcount[,2]
alt_SH <- pooldata@readcoverage[,3] - pooldata@refallele.readcount[,3]
alt_ARA <- pooldata@readcoverage[,4] - pooldata@refallele.readcount[,4]
alt_BMR <- pooldata@readcoverage[,5] - pooldata@refallele.readcount[,5]
alt_CBL <- pooldata@readcoverage[,6] - pooldata@refallele.readcount[,6]
alt_FC <- pooldata@readcoverage[,7] - pooldata@refallele.readcount[,7]
alt_FR <- pooldata@readcoverage[,8] - pooldata@refallele.readcount[,8]
alt_HZD <- pooldata@readcoverage[,9] - pooldata@refallele.readcount[,9]
alt_SBR <- pooldata@readcoverage[,10] - pooldata@refallele.readcount[,10]
alt_PSG <- pooldata@readcoverage[,11] - pooldata@refallele.readcount[,11]
alt_KH <- pooldata@readcoverage[,12] - pooldata@refallele.readcount[,12]
alt_STC <- pooldata@readcoverage[,13] - pooldata@refallele.readcount[,13]
alt_PL <- pooldata@readcoverage[,14] - pooldata@refallele.readcount[,14]
alt_VD <- pooldata@readcoverage[,15] - pooldata@refallele.readcount[,15]
alt_OCT <- pooldata@readcoverage[,16] - pooldata@refallele.readcount[,16]
alt_PB <- pooldata@readcoverage[,17] - pooldata@refallele.readcount[,17]
alt_PGP <- pooldata@readcoverage[,18] - pooldata@refallele.readcount[,18]
alt_PSN <- pooldata@readcoverage[,19] - pooldata@refallele.readcount[,19]

# Calculate frequency for each population
fq_STR <- ref_STR/pooldata@readcoverage[,1]
fq_SLR <- ref_SLR/pooldata@readcoverage[,2]
fq_SH <- ref_SH/pooldata@readcoverage[,3]
fq_ARA <- ref_ARA/pooldata@readcoverage[,4]
fq_BMR <- ref_BMR/pooldata@readcoverage[,5]
fq_CBL <- ref_CBL/pooldata@readcoverage[,6]
fq_FC <- ref_FC/pooldata@readcoverage[,7]
fq_FR <- ref_FR/pooldata@readcoverage[,8]
fq_HZD <- ref_HZD/pooldata@readcoverage[,9]
fq_SBR <- ref_SBR/pooldata@readcoverage[,10]
fq_PSG <- ref_PSG/pooldata@readcoverage[,11]
fq_KH <- ref_KH/pooldata@readcoverage[,12]
fq_STC <- ref_STC/pooldata@readcoverage[,13]
fq_PL <- ref_PL/pooldata@readcoverage[,14]
fq_VD <- ref_VD/pooldata@readcoverage[,15]
fq_OCT <- ref_OCT/pooldata@readcoverage[,16]
fq_PB <- ref_PB/pooldata@readcoverage[,17]
fq_PGP <- ref_PGP/pooldata@readcoverage[,18]
fq_PSN <- ref_PSN/pooldata@readcoverage[,19]


# pooldata@snp.info is a data frame (nsnp row and 4 columns) detailing for each SNP, the chromosome (or scaffold), the position, Reference allele name and Alternate allele name

# Name each SNP based on its chromosome/scaffold and position
SNP <- paste(pooldata@snp.info[,1], pooldata@snp.info[,2], sep="_")

# Create matrix and rename columns and rows
poolfstat_matrix = matrix(nrow=19, ncol=length(SNP))
colnames(poolfstat_matrix) <- SNP
rownames(poolfstat_matrix)=c("STR","SLR","SH","ARA", "BMR", "CBL", "FC", "FR", "HZD", "SBR", "PSG", "KH", "STC", "PL", "VD", "OCT", "PB", "PGP", "PSN")

# Fill matrix with relative frequencies for each population
poolfstat_matrix[1,]=fq_STR
poolfstat_matrix[2,]=fq_SLR
poolfstat_matrix[3,]=fq_SH
poolfstat_matrix[4,]=fq_ARA
poolfstat_matrix[5,]=fq_BMR
poolfstat_matrix[6,]=fq_CBL
poolfstat_matrix[7,]=fq_FC
poolfstat_matrix[8,]=fq_FR
poolfstat_matrix[9,]=fq_HZD
poolfstat_matrix[10,]=fq_SBR
poolfstat_matrix[11,]=fq_PSG
poolfstat_matrix[12,]=fq_KH
poolfstat_matrix[13,]=fq_STC
poolfstat_matrix[14,]=fq_PL
poolfstat_matrix[15,]=fq_VD
poolfstat_matrix[16,]=fq_OCT
poolfstat_matrix[17,]=fq_PB
poolfstat_matrix[18,]=fq_PGP
poolfstat_matrix[19,]=fq_PSN

# ================================================================================== #
# ================================================================================== #

# Run pcadapt

# Load matrix into pcadapat
pool.data <- read.pcadapt(poolfstat_matrix, type = "pool")

# Read in metadata
meta_path <- "data/processed/outlier_analyses/guide_files/Populations_metadata.csv"
meta <- read.csv(meta_path, header=T)

# ================================================================================== #

# With Pool-Seq data, the package computes again a Mahalanobis distance based on PCA loadings.

# Note: To choose K, principal component analysis should first be performed with a large enough number of principal components. 
# Use scree plot and PCA to identify optimal K. 

# Run pcadapt on data with large number of K
# You can also set the parameter min.maf that corresponds to a threshold of minor allele frequency. 
# By default, the parameter min.maf is set to 5%
x <- pcadapt(input = pool.data, K = 18)

# Graph screeplot (The ideal pattern in a scree plot is a steep curve followed by a bend and a straight line. 
# The eigenvalues that correspond to random variation lie on a straight line whereas the ones that correspond to population structure 
# lie on a steep curve). Keep PCs that correspond to eigenvalues to the left of the straight line (Cattell’s rule)
pdf("output/figures/outlier_analyses/pcadapt_screeplot_K18.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "screeplot")
dev.off()

# Graph PCA - PC1 and PC2
pdf("output/figures/outlier_analyses/pcadapt_pca_PC1_PC2_K18.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "scores")
dev.off()

# Graph PCA - PC3 and PC4
pdf("output/figures/outlier_analyses/pcadapt_pca_PC3_PC4_K18.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "scores", i = 3, j = 4)
dev.off()

# ================================================================================== #

# Run pcadapt on data with optimal K
x <- pcadapt(input = pool.data, K = 5)
pc.percent<-round(100*(x$singular.values^2), digits = 2)
# 69.01  9.67  6.89  3.42  2.95

# Summary of pcadapt
summary(x)

# ================================================================================== #

# Identify outliers 

# Scree plot
pdf("output/figures/outlier_analyses/pcadapt_screeplot_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "screeplot")
dev.off()

# Color palette 
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))
colors.reorder <- mycolors[c(19,2,3,4,11,5,1,10,17,14,6,8,7,13,9,18,16,12,15)]
colors.alphabetical <- mycolors[c(4,11,5,1,10,17,8,18,16,12,13,6,15,14,3,2,7,19,9)]

# Graph PCA - PC1 and PC2
pdf("output/figures/outlier_analyses/pcadapt_pca_PC1_PC2_K5.pdf", width = 9, height = 8)
pcadapt:::plot.pcadapt(x, option = "scores", pop=rownames(poolfstat_matrix), col=colors.alphabetical) 
dev.off()
# Note: need to specify pop. If field is left empty, the points will be displayed in black.

# Graph PCA - PC3 and PC4
pdf("output/figures/outlier_analyses/pcadapt_pca_PC3_PC4_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "scores", i = 3, j = 4, pop=rownames(poolfstat_matrix), col=colors.alphabetical)
dev.off()

# Manhattan plot
pdf("output/figures/outlier_analyses/pcadapt_manhattan_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "manhattan")
dev.off()

# QQ plot
pdf("output/figures/outlier_analyses/pcadapt_qqplot_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x, option = "qqplot")
dev.off()

# Histogram of test statistic and p-val 
pdf("output/figures/outlier_analyses/pcadapt_hist_K5.pdf", width = 8, height = 8)
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange")
dev.off()

# Histogram of the test statistic 𝐷𝑗
pdf("output/figures/outlier_analyses/pcadapt_stat.dist_K5.pdf", width = 8, height = 8)
plot(x, option = "stat.distribution")
dev.off()

# ================================================================================== #

# Choose cut-off for outlier detection

# q-value method 
# SNPs with q-values less than 𝛼 (10%) will be considered as outliers with an expected false discovery rate bounded by 𝛼
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.001
outliers <- which(qval < alpha)
length(outliers) #794,830
snp_outliers_pc <- get.pc(x, outliers) # Get outlier snps associated with PCs
outliers <- data.frame(outliers)
outliers$line_num <- 1:nrow(outliers) # Add line numbers

# Benjamini-Hochberg Procedue
padj_BH <- p.adjust(x$pvalues, method="BH")
alpha <- 0.001
outliers_BH <- which(padj_BH < alpha)
length(outliers_BH) #794,830
snp_outliers_BH_pc <- get.pc(x, outliers_BH) # Get outlier snps associated with PCs
outliers_BH <- data.frame(outliers_BH)
outliers_BH$line_num <- 1:nrow(outliers_BH) # Add line numbers

# Bonferroni Correction
padj_bonferroni <- p.adjust(x$pvalues, method="bonferroni")
alpha <- 0.001
outliers_bonferroni <- which(padj_bonferroni < alpha)
length(outliers_bonferroni) #25,103
snp_outliers_bonferroni_pc <- get.pc(x, outliers_bonferroni) # Get outlier snps associated with PCs
outliers_bonferroni <- data.frame(outliers_bonferroni) 
outliers_bonferroni$line_num <- 1:nrow(outliers_bonferroni) # Add line numbers

# ================================================================================== #

# Save outlier SNP lists
write.table(outliers, file = "data/processed/outlier_analyses/pcadapt_outliers_results.snp", quote = F, row.names = F, col.names = F)
write.table(outliers_BH, file = "data/processed/outlier_analyses/pcadapt_outliers_BH_results.snp", quote = F, row.names = F, col.names = F)
write.table(outliers_bonferroni, file = "data/processed/outlier_analyses/pcadapt_outliers_bonferroni_results.snp", quote = F, row.names = F, col.names = F)

# ================================================================================== #

# Try thinning to remove LD
# size = window radius 
# thr = the squared correlation threshold
x_LD <- pcadapt(input = pool.data, K = 5, LD.clumping = list(size = 500, thr = 0.1))

# Summary
summary(x_LD)

# ================================================================================== #

# Scree plot
pdf("output/figures/outlier_analyses/pcadapt_screeplot_LD_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x_LD, option = "screeplot")
dev.off()

# Graph PCA - PC1 and PC2
pdf("output/figures/outlier_analyses/pcadapt_pca_PC1_PC2_LD_K5.pdf", width = 9, height = 8)
pcadapt:::plot.pcadapt(x_LD, option = "scores", pop=rownames(poolfstat_matrix), col=colors.alphabetical) 
dev.off()

# Graph PCA - PC3 and PC4
pdf("output/figures/outlier_analyses/pcadapt_pca_PC3_PC4_LD_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x_LD, option = "scores", i = 3, j = 4, pop=rownames(poolfstat_matrix), col=colors.alphabetical)
dev.off()

# Display the loadings (contributions of each SNP to the PC) to evalutate if the loadings are clustered in a single or several genomic regions
# You want the laodings to be evenly distributed 
pdf("output/figures/outlier_analyses/pcadapt_loadings1:4.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
for (i in 1:4)
  plot(x_LD$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))
dev.off()

# Manhattan plot
pdf("output/figures/outlier_analyses/pcadapt_manhattan_LD_K5.pdf", width = 8, height = 8)
pcadapt:::plot.pcadapt(x_LD, option = "manhattan")
dev.off()