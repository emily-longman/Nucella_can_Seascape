# Fisher's Exact Test of SNP types

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
#install.packages(c('data.table', 'tidyverse', 'foreach', 'dplyr', 'RColorBrewer'))
library(data.table)
library(tidyverse)
library(foreach)
library(dplyr)
library(RColorBrewer)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("output/figures/outlier_analyses")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Read in SNP data
snp.meta <- read.table("data/processed/baypass/input_files/snpdet", header=F)
# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")

# Load data - outliers for pH and Mcali based on beating POD BF threshold
bf.ph.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header = T)
bf.McaliIntThk.mean.sum.outliers.annotated <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", header = T)

# Load annotation
load("data/processed/outlier_analyses/Ncan.pooldata.annotations.RData")

# ================================================================================== #

# Summarize each list based on annotation
ph_ann_sum <- bf.ph.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(ph = n)
McaliIntThk_ann_sum <- bf.McaliIntThk.mean.sum.outliers.annotated %>% count(Annotation) %>% rename(McaliIntThk = n)

# Join data
ann <- full_join(ph_ann_sum, McaliIntThk_ann_sum)

# Summarize full SNP list
all_ann <- o.merge.filt %>% count(col) %>% rename(Annotation = col, genome = n)

# Join data
ann <- full_join(ann, all_ann)
ann[is.na(ann)] <- 0

#--------------------------------------------------------------------------------

# Subset data for common var
ann.focal <- c("missense_variant", "synonymous_variant", "intron_variant", "intergenic_region", 
"upstream_gene_variant", "downstream_gene_variant", "3_prime_UTR_variant", "5_prime_UTR_variant")
ann.sub <- ann %>% filter(Annotation %in% ann.focal)

#--------------------------------------------------------------------------------

# Loop through annotations and perform Fishers exact test - pH
ftests_ph <- foreach(i=ann.focal, .combine = "rbind", .errorhandling = "remove")%do%{
  
  # Create matrix for significance of focal annotation
  ann.tmp <- matrix(c(
            ann$ph[which(ann$Annotation==i)],
            sum(ann$ph)-ann$ph[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$ph[which(ann$Annotation==i)],
            sum(ann$genome)-ann$genome[which(ann$Annotation==i)]),
            nrow = 2)

  # Fishers exact test
  ftest_i <- fisher.test(ann.tmp)

  # Create data frame with output (i, p.value, odds ratio, and 95% CI)
  data.frame(
      class = i,
      p.fet = ftest_i$p.value,
      OR = ftest_i$estimate,
      lci = ftest_i$conf.int[1],
      uci = ftest_i$conf.int[2]
    )
}

# Prep for graphing
ftests_ph$class_graphing <- ftests_ph$class
ftests_ph$class_graphing <- str_remove_all(ftests_ph$class_graphing, "_variant")
ftests_ph$class_graphing <- gsub("_", " ", ftests_ph$class_graphing)

# Order classes based on OR
ftests_ph <- ftests_ph %>% mutate(class_graphing = fct_reorder(class_graphing, OR))

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test_pH.pdf", width = 12, height = 8)
ggplot(ftests_ph, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet))) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1) +
  geom_point(shape = 21, size = 8) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") +
  theme_bw(base_size=36) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 22), legend.position = c(0.86, 0.22), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()


# Loop through annotations and perform Fishers exact test - M californianus
ftests_Mcali <- foreach(i=ann.focal, .combine = "rbind", .errorhandling = "remove")%do%{
  
  # Create matrix for significance of focal annotation
  ann.tmp <- matrix(c(
            ann$McaliIntThk[which(ann$Annotation==i)],
            sum(ann$McaliIntThk)-ann$McaliIntThk[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$McaliIntThk[which(ann$Annotation==i)],
            sum(ann$genome)-ann$genome[which(ann$Annotation==i)]),
            nrow = 2)

  # Fishers exact test
  ftest_i <- fisher.test(ann.tmp)

  # Create data frame with output (i, p.value, odds ratio, and 95% CI)
  data.frame(
      class = i,
      p.fet = ftest_i$p.value,
      OR = ftest_i$estimate,
      lci = ftest_i$conf.int[1],
      uci = ftest_i$conf.int[2]
    )
}

# Prep for graphing
ftests_Mcali$class_graphing <- ftests_Mcali$class
ftests_Mcali$class_graphing <- str_remove_all(ftests_Mcali$class_graphing, "_variant")
ftests_Mcali$class_graphing <- gsub("_", " ", ftests_Mcali$class_graphing)

# Order classes based on OR
ftests_Mcali <- ftests_Mcali %>% mutate(class_graphing = fct_reorder(class_graphing, OR))

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test_Mcali.pdf", width = 12, height = 8)
ggplot(ftests_Mcali, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet))) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1) +
  geom_point(shape = 21, size = 8) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") + 
  #labs(title = expression(paste(italic("M. californianus"), " cross-sectional thickness"))) +
  theme_bw(base_size=36) + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10,, unit = "pt")) +
  theme(legend.title = element_text(size = 24), legend.text = element_text(size = 22), legend.position = c(0.15, 0.81), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

#--------------------------------------------------------------------------------

# Graph both pH and Mcali in same figure

# Merge for graphing
ftests_ph$model <- "ph"
ftests_Mcali$model <- "Mcali"

# Join
ftests_all <- rbind(ftests_ph, ftests_Mcali)

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test_all.pdf", width = 10, height = 6)
ggplot(ftests_all, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet), shape=model)) + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_point(size = 5, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") +
  #scale_fill_gradientn(colours=brewer.pal(9, "Greens"), name="-log10(p)") + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") +
  theme_bw(base_size=22) 
dev.off()


#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# Directly compare pH and Mcali models

# Loop through annotations and perform Fishers exact test
ftests <- foreach(i=ann.focal, .combine = "rbind", .errorhandling = "remove")%do%{
  
  # Create matrix for significance of focal annotation
  ann.tmp <- matrix(c(
            ann$ph[which(ann$Annotation==i)],
            ann$McaliIntThk[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$ph[which(ann$Annotation==i)],
            ann$genome[which(ann$Annotation==i)]-ann$McaliIntThk[which(ann$Annotation==i)]),
            nrow = 2)

  # Fishers exact test
  ftest_i <- fisher.test(ann.tmp)

  # Create data frame with output (i, p.value, odds ratio, and 95% CI)
  data.frame(
      class = i,
      p.fet = ftest_i$p.value,
      OR = ftest_i$estimate,
      lci = ftest_i$conf.int[1],
      uci = ftest_i$conf.int[2]
    )
}

# Create new column for graphing
ftests$class_graphing <- ftests$class
ftests$class_graphing <- str_remove_all(ftests$class_graphing, "_variant")
ftests$class_graphing <- gsub("_", " ", ftests$class_graphing)

# Order classes based on OR
ftests <- ftests %>% mutate(class_graphing = fct_reorder(class_graphing, OR))

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test.pdf", width = 10, height = 6)
ggplot(ftests, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet))) + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1) +
  geom_point(shape = 21, size = 5) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") +
  #scale_fill_gradientn(colours=brewer.pal(9, "Greens"), name="-log10(p)") + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") +
  theme_bw(base_size=22)
dev.off()
