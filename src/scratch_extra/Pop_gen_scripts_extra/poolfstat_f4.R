# Calculate F4

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
install.packages(c('poolfstat', 'tidyverse', 'ggplot2', 'RColorBrewer', 'maps', 'mapdata', 'ggrepel', 'pheatmap'))
library(poolfstat)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(maps) 
library(mapdata)
library(ggrepel)
library(pheatmap)

# ================================================================================== #

# Read in population names
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)

# ================================================================================== #

# Create a pooldata object for Pool-Seq read count data (poolsize = haploid sizes of each pool, # of pools)
# Note: 20 individuals per pool. N. canaliculata is a diploid species. So haploid size = 40 for most pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 26,206,958)
#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 11,656,080 SNPs for 19 Pools

# Read in data and filter (note: the input vcf is the vcf output from 12_filter_vcf.sh - # variants: 14,897,468)
pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf", 
poolsizes=rep(40,19), poolnames=pops$V1, 
min.cov.per.pool = 20, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 8,277,206 SNPs for 19 Pools

#pooldata <-vcf2pooldata(vcf.file="data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ50_maxmissing0.9.recode.vcf", 
#poolsizes=rep(40,19), poolnames=pops$V1, 
#min.cov.per.pool = 15, min.rc = 5, max.cov.per.pool = 120, min.maf = 0.01, nlines.per.readblock = 1e+06)
# Data consists of 9,804,957 SNPs for 19 Pools

# min.cov.per.pool = the minimum allowed read count per pool for SNP to be called
# min.rc =  the minimum # reads that an allele needs to have (across all pools) to be called 
# max.cov.per.pool = the maximum read count per pool for SNP to be called 
# min.maf = the minimum allele frequency (over all pools) for a SNP to be called (note this is obtained from dividing the read counts for the minor allele over the total read coverage) 
# nlines.per.readblock = number of lines in sync file to be read simultaneously 

# Save pooldata
save(pooldata, file="data/processed/pop_structure/pooldata.RData")
# Reload pooldata
load("data/processed/pop_structure/pooldata.RData")

# ================================================================================== #


# Compute parameters Fs, F3, F4, D, heterozygosities

# Estimation of f-statistics on Pool-Seq data (with computation of Dstat)
pooldata.fstats <- compute.fstats(pooldata, nsnp.per.bjack.block = 1000, return.F2.blockjackknife.samples = TRUE)

# Save pooldata
save(pooldata.fstats, file="data/processed/pop_structure/pooldata.fstats.RData")
# Reload pooldata.fstats
load("data/processed/pop_structure/pooldata.fstats.RData")

# ================================================================================== #

# Test for four-population treeness (note:  a Z-score lower than 1.96 in absolute value provides no evidence against the null-hypothesis of treeness for the tested population configuration at the 95% significance threshold)
tst.sel<-abs(pooldata.fstats@f4.values$`Z-score`)<1.96
as.data.frame(pooldata.fstats@f4.values)[tst.sel,]

# Estimating admixture proportions with f4-ratios
# Example of estimating f4 ratios
compute.f4ratio(pooldata.fstats, num.quadruplet = "HZD,PSN;PB,OCT", den.quadruplet="HZD,PSN;SBR,STR")
#  Estimate bjack mean bjack s.e.    CI95inf    CI95sup 
# 3.0296472  2.9891599  0.0992849  2.7945615  3.1837583 
# Note: Simulated value (α = 0.25) is within the confidence interval

# ================================================================================== #

# Build admixture graph from scratch

# The find.tree.popset function selects maximal sets of unadmixed populations from an fstats object
scaf.pops <- find.tree.popset(pooldata.fstats, verbose=FALSE)

# List the 15 passing the treeness test for the identified set
scaf.pops$passing.quaduplet

# State the range of variation of the passing quadruplet
scaf.pops$Z_f4.range

# The rooted.njtree.builder function builds the scaffold tree based on the set of unadmixed populations (identified above with find.tree.popset) 
scaf.tree <- rooted.njtree.builder(fstats=pooldata.fstats, pop.sel=scaf.pops$pop.sets[1,], plot.nj=FALSE)
# Score of the NJ tree: 27.19838 

# Check the fit of the neighbor-joining tree
scaf.tree$nj.tree.eval

# Note: need to graph using R GUI

# Plot best inferred rooted scaffold tree
plot(scaf.tree$best.rooted.tree)
# This tree can be used as a reference graph to construct the complete admixture graph

# Add other populations
add.pool.OCT <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="OCT", 
                         fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT$best.fitted.graph)

add.pool.OCT.HZD <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="HZD", 
                         fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT.HZD$best.fitted.graph)

add.pool.OCT.HZD.SBR <- add.leaf(scaf.tree$best.rooted.tree,leaf.to.add="SBR", 
                             fstats=pooldata.fstats, verbose=FALSE, drift.scaling=TRUE)

plot(add.pool.OCT.HZD.SBR$best.fitted.graph)