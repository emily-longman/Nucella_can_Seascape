# Use poolfstat to convert VCF to Baypass input file

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
install.packages(c('poolfstat', 'tidyverse', 'ggplot2', 'RColorBrewer', 'viridis'))
library(poolfstat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(viridis)

# ================================================================================== #

# Generate Folders and files

# Make data directory
data_dir="data/raw/pooldata"
if (!dir.exists(data_dir)) {dir.create(data_dir)}

# Make output directory
output_dir="output/figures/demography"
if (!dir.exists(output_dir)) {dir.create(output_dir)}

# ================================================================================== #

# Read in population names
pops <- read.table("data/raw/vcf_clean/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter
pooldata <-vcf2pooldata(vcf.file="data/raw/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# ================================================================================== #

# Save pooldata
save(pooldata, file="data/raw/pooldata/pooldata.RData")
# Reload pooldata
load("data/raw/pooldata/pooldata.RData")

# ================================================================================== #

# Principle Components Analysis with randomallele.pca

# PCA on the read count data (the object)
pooldata.pca = randomallele.pca(pooldata, main="Read Count data")

# Save loadings as data frame
pooldata.pca$pop.loadings %>% as.data.frame -> pca.df

# Rename columns
colnames(pca.df) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18")

# ================================================================================== #

# Save dataframe
write.csv(pca.df, "data/processed/outlier_analyses/pca.csv")

# ================================================================================== #

# Graph PCA

# Add column that highlights N, and S
pca.df <- pca.df %>% mutate(site = rownames(pca.df), 
    shape = case_when(site %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   site %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))

# Add column that highlights N, S, admixture
#pca.df <- pca.df %>% mutate(site = rownames(pca.df), 
#    shape = case_when(site %in% c("STR", "OCT", "HZD", "PB", "PSN") ~ "S", 
#                   site %in% c("SBR", "PL", "PGP") ~ "Admix", 
#                   site %in% c("BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))

# Order by N to S
lat <- c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR")
pca.df.order <- pca.df %>% mutate(site = factor(site, levels = lat)) %>% arrange(site)

# Color palette
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(19))
#colors.reorder <- mycolors[c(19,2,3,4,11,5,1,10,17,14,6,8,7,13,9,18,16,12,15)]
viridiscolors <- viridis(n=19)

# Plot PC1 and PC2 with ggplot
pdf("output/figures/demography/PCA_all_SNPs_PC1_PC2_ggplot_bigger_NvsS.pdf", width = 12, height = 9)
ggplot(pca.df.order, aes(x=PC1, y=PC2, shape=shape, fill = factor(site))) + geom_jitter(size=16, width = 0.015, height = 0.015) + 
scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
ylim(-0.51, 0.51) + xlim(-0.51, 0.51) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) + xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) +
theme_linedraw(base_size = 36) +
geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none") +
labs(fill = "Site")
dev.off()

# Plot with base R
pdf("output/figures/demography/PCA_all_SNPs_PC1_PC2.pdf", width = 8, height = 8)
par(mar=c(5,6,4,1)+.1) # Adjust margins
pca <- plot(pooldata.pca$pop.loadings[,1],pooldata.pca$pop.loadings[,2],
xlab=paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)"),
ylab=paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)"),
col="black", bg=colors.reorder, pch=21, cex = 5, cex.lab = 3)
abline(h=0,lty=2,col="grey"); abline(v=0,lty=2,col="grey")
dev.off()

