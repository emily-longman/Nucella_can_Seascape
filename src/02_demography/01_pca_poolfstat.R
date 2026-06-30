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
ggplot(pca.df.order, aes(x=PC1, y=PC2, shape=shape, fill = factor(site))) + geom_jitter(size=16, width = 0.01, height = 0.01) + 
scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
ylim(-0.51, 0.51) + xlim(-0.51, 0.51) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) + xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) +
theme_linedraw(base_size = 40) +
geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none") +
labs(fill = "Site")
dev.off()
pdf("output/figures/demography/PCA_all_SNPs_PC1_PC2_ggplot_bigger_NvsS_alt.pdf", width = 12, height = 9)
ggplot(pca.df.order, aes(x=PC1, y=PC2, shape=shape, fill = factor(site))) + geom_jitter(size=16, width = 0.01, height = 0.01, alpha=0.8) + 
scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
ylim(-0.51, 0.51) + xlim(-0.51, 0.51) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) + xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) +
theme_linedraw(base_size = 40) +
geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none") +
labs(fill = "Site")
dev.off()
pdf("output/figures/demography/PCA_all_SNPs_PC1_PC2_ggplot_bigger_NvsS_nojitter.pdf", width = 12, height = 9)
ggplot(pca.df.order, aes(x=PC1, y=PC2, shape=shape, fill = factor(site))) + geom_point(size=16) + 
scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
ylim(-0.51, 0.51) + xlim(-0.51, 0.51) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) + xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) +
theme_linedraw(base_size = 40) +
geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none") +
labs(fill = "Site")
dev.off()
pdf("output/figures/demography/PCA_all_SNPs_PC1_PC2_ggplot_bigger_NvsS_nojitter_alt.pdf", width = 12, height = 9)
ggplot(pca.df.order, aes(x=PC1, y=PC2, shape=shape, fill = factor(site))) + geom_point(size=16, alpha=0.7) + 
scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
ylim(-0.51, 0.51) + xlim(-0.51, 0.51) + 
ylab(paste0("PC",2," (",round(pooldata.pca$perc.var[2],2),"%)")) + xlab(paste0("PC",1," (",round(pooldata.pca$perc.var[1],2),"%)")) +
theme_linedraw(base_size = 40) +
geom_vline(xintercept = 0, color = "black", linetype = "dashed") + geom_hline(yintercept = 0, color = "black", linetype = "dashed") + 
guides(fill = guide_legend(override.aes = list(shape = c(21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 21, 23, 23, 23, 23, 23, 23, 23), size = 8)), shape = "none") +
labs(fill = "Site")
dev.off()

# ================================================================================== #

# Graph map of sites

# Get state data
states <- map_data("state")
# Subset data for only California and Oregon
west_coast <- subset(states, region %in% c("california", "oregon"))

# Coordinates of sites
sites <- data.frame(
  longitude = c(-124.0593, -124.0848, -124.1148, -124.4015, -124.5647, -124.2529, -124.0809, -123.7895, -123.8036, 
                -123.2551, -123.0740, -122.3976, -121.9537, -121.9290, -121.3187, -121.2868, -120.8838, -120.6399, -120.6157),
  latitude = c(44.83777, 44.50540, 44.24999, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461, 39.28090, 38.51198, 38.31900, 
               37.18506, 36.51939, 36.44750, 35.72893, 35.66549, 35.28994, 34.88117, 34.73024),
  site.abrev = c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

sites <- sites %>% mutate(shape = case_when(site.abrev %in% c("STR", "OCT", "HZD", "PB", "PSN", "SBR", "PL") ~ "S", 
                   site.abrev %in% c("PGP", "BMR", "FR", "VD", "KH", "STC", "PSG", "CBL", "ARA", "SH", "SLR", "FC") ~ "N"))
sites <- sites %>% mutate(site.abrev = factor(site.abrev, levels = lat)) %>% arrange(site.abrev)


# Coordinates of site labels
lat.site.labels <- c(44.83777+0.09, 44.50540+0.02, 44.24999-0.09, 43.30402, 42.84097, 41.77121, 40.03011, 39.60461-0.04, 39.28090-0.07, 38.51198+0.05, 38.31900-0.11, 
               37.18506+0.05, 36.51939+0.15, 36.44750-0.1, 35.72893+0.17, 35.66549-0.1, 35.28994-0.07, 34.88117-0.08, 34.73024-0.27)

long.site.labels.abrev <- c(-124.0593-0.59, -124.0848-0.75, -124.1148-0.59, -124.4015-0.75, -124.5647-0.72, -124.2529-0.75, 
                      -124.0809-0.75, -123.7895-0.5, -123.8036-0.64, -123.2551-0.6, -123.0740-0.75, -122.3976-0.75, 
                      -121.9537-0.7, -121.9290-0.75, -121.3187-0.78, -121.2868-0.75, -120.8838-0.72, -120.6399-0.78, -120.6157-0.65)


# Create site map (site codes)
pdf("output/figures/Site_map.pdf", width = 8, height = 9)
ggplot(data = west_coast) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude, shape = shape, fill = factor(site.abrev)), size = 8) + coord_fixed(1.3) +
  scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
  geom_text(data=sites, aes(long.site.labels.abrev-0.4, lat.site.labels, label=site.abrev), size = 6)+ 
  scale_x_continuous(limits = c(-126, -114), breaks = seq(-125, -114, by = 3)) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + theme(legend.position = "none", panel.grid = element_blank())
dev.off()

pdf("output/figures/Site_map.pdf", width = 8, height = 14)
ggplot(data = west_coast) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude, shape = shape, fill = factor(site.abrev)), size = 10) + coord_fixed(1.3) +
  scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = viridiscolors) +
  geom_text(data=sites, aes(long.site.labels.abrev-0.4, lat.site.labels, label=site.abrev), size = 8)+ 
  scale_x_continuous(limits = c(-126, -114.1), breaks = seq(-125, -114.1, by = 3), expand = expansion(mult = c(0.05, 0.01))) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 30) + theme(legend.position = "none", panel.grid = element_blank())
dev.off()

# Join sites with PCs
pca.df.tmp <- pca.df %>% rename(site.abrev = site)
sites <- left_join(sites, pca.df.tmp)

pdf("output/figures/Site_map_PC1.pdf", width = 9.5, height = 14)
ggplot(data = west_coast) +
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "black") + 
  geom_point(data = sites, aes(x = longitude, y = latitude, shape = shape, fill = PC1), size = 11) + coord_fixed(1.3) +
  scale_shape_manual(values = c(21, 23)) + #scale_fill_manual(values = viridiscolors) +
  scale_fill_gradientn(colours=brewer.pal(9, "RdGy"), breaks = c(-0.2, 0.00)) +
  geom_text(data=sites, aes(long.site.labels.abrev-0.35, lat.site.labels, label=site.abrev), size = 8) + 
  scale_x_continuous(limits = c(-125.75, -114.1), breaks = seq(-125, -114.1, by = 3), expand = expansion(mult = c(0.04, 0.01))) +
  xlab("Longitude") + ylab("Latitude") + theme_linedraw(base_size = 32) + 
  theme(panel.grid = element_blank(), legend.title = element_text(size = 32), legend.text = element_text(size = 30), legend.key.size = unit(1.2, "cm"), legend.position = c(0.81, 0.53)) + 
  guides(shape = "none")
dev.off()
