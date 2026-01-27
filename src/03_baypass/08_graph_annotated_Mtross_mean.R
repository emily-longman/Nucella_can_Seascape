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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'ggplot2', 'RColorBrewer', 'poolfstat))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(poolfstat)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/baypass")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

# Load baypass BF data (bf.Mtross.mean.sum)
load("data/processed/baypass/biotic/bf.Mtross.mean.sum.Rdata")

# Load POD thresholds
load("data/processed/baypass/biotic/Mtross_mean_POD_thr.Rdata")

# Load annotated outliers
bf.Mtross.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.Mtross.mean.sum.outliers.annotated.csv", header=T)

# ================================================================================== #

# How many unique genes - 2237
length(unique(bf.Mtross.mean.sum.outliers.annotated$Gene_Name))

# Identify points in gene g26813
g26813 <- bf.Mtross.mean.sum.outliers.annotated %>% filter(Gene_Name == "g26813")
Mtross <- bf.Mtross.mean.sum %>% mutate(g26813 = MRK %in% g26813$MRK)


# Graph BF with 0.001 POD threshold - only BF > 0 and coloring points in gene g26813
pdf("output/figures/baypass/baypass_BF_Mtross_mean_repmeans_posBF_g26813.pdf", width = 12, height = 8)
ggplot(Mtross[which(Mtross$bf_db.mean>0),], aes(y=bf_db.mean, x=chr, color=g26813)) + 
  labs(x = "Position", y = "BF (in dB)") +
  scale_color_manual(values=c("#bebebecc", "black")) +
  geom_point(alpha=0.6) + 
  geom_hline(yintercept=bf.POD.thr$bf_db.mean[which(bf.POD.thr$thr==0.999)], linetype = "dashed", col="red") +
  theme_classic(base_size = 20) +
  theme(legend.position = "none") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 14))
dev.off()

# ================================================================================== #

# Graph against latitude

# Load pooldata
load("data/raw/pooldata/pooldata.RData")

# Get the snp.info matrix
snp.info <- pooldata@snp.info
snp.info <- snp.info %>% mutate(SNP_id = paste(Chromosome, Position, sep = "_"))
# Add column for g26813
snp.info <- snp.info %>% mutate(g26813 = SNP_id %in% g26813$SNP_id)


# Subset the pooldata object using the selected SNP indices
pooldata.sub <- pooldata.subset(pooldata, snp.index = which(snp.info$g26813 == TRUE))



# Extract read count and coverage data for SNPs
snp.info.sub <- pooldata.sub@snp.info
ref_count <- pooldata.sub@refallele.readcount
coverage <- pooldata.sub@readcoverage

# Calculate allele frequency for SNPs
afs <- ref_count/coverage
afs <- afs %>% as.data.frame()

# Rename columns (19 sites and SNP_ids)
names(afs) = c(pooldata.sub@poolnames)
rownames(afs) = rownames(snp.info.sub)

# Graph allele frequencies

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# AF versus latitude
pdf("output/figures/baypass/Outlier_AF_lat.pdf", width = 18, height = 4)
ggplot(data = afs.id.mapped.targets, aes(x = Lat, y = AF)) +
  geom_point(aes(fill = Site), size = 5, shape = 21) + 
  geom_vline(xintercept=36.8007, linetype="solid", color="black") +
  ylab("Allele Frequency") + xlab("Latitude") +
  scale_fill_manual(values = mycolors) +
  facet_grid(~pos) + theme_bw(base_size = 22) + theme(legend.position="none")
dev.off()


# Graph map
pdf("output/figures/morphology/Allele_freq_maps.pdf", width = 8, height = 3)
ggplot(data = world) +
  geom_sf(fill= "grey70") +
  coord_sf(xlim = c(-127, -115), ylim = c(32, 47), expand = FALSE) + 
  theme(panel.grid.major = element_line(color = gray(.5), linetype = "dashed", linewidth = 0.2), 
        panel.background = element_rect(fill = "white")) +
  geom_jitter(data = afs.id.mapped.targets,
              color = "black",
              aes(
                x=Long,
                y=Lat,
                fill = AF,
                shape = cluster
              ), alpha = 0.9, size = 2.5) + 
  scale_fill_gradient2(low = "steelblue", ,high = "firebrick",
                       midpoint = 0.5) +
  scale_shape_manual(values = 21:23) +
  facet_grid(~snp_id) + theme_bw()
dev.off()
