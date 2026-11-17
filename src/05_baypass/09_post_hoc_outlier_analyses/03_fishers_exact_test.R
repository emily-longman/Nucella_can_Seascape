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
bf.ph.mean.sum.outliers.annotated_multiann <- read.csv("data/processed/baypass/bf.ph.mean.sum.outliers.annotated.csv", header = T)
bf.McaliIntThk.mean.sum.outliers.annotated_multiann <- read.csv("data/processed/baypass/bf.McaliIntThk.mean.sum.outliers.annotated.csv", header = T)

# Load annotation
load("data/processed/outlier_analyses/Ncan.pooldata.annotations.RData")

# ================================================================================== #

# Keep only top annotation
bf.ph.mean.sum.outliers.annotated <- bf.ph.mean.sum.outliers.annotated_multiann %>% filter(annotation.id == 1)
bf.McaliIntThk.mean.sum.outliers.annotated <- bf.McaliIntThk.mean.sum.outliers.annotated_multiann  %>% filter(annotation.id == 1)

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

# Save table
write.csv(ftests_ph, "output/tables/ftests_ph.csv", row.names=F)

# Prep for graphing
ftests_ph$class_graphing <- ftests_ph$class
ftests_ph$class_graphing <- str_remove_all(ftests_ph$class_graphing, "_variant")
ftests_ph$class_graphing <- gsub("_", " ", ftests_ph$class_graphing)

# Order classes based on OR
ftests_ph <- ftests_ph %>% mutate(class_graphing = fct_reorder(class_graphing, OR))

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test_pH.pdf", width = 8, height = 8)
ggplot(ftests_ph, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet))) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1) +
  geom_point(shape = 21, size = 8) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") +
  theme_linedraw(base_size=30) + theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.title = element_text(size = 19), legend.text = element_text(size = 20), legend.position = c(0.8, 0.178), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

#--------------------------------------------------------------------------------

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

# Save table
write.csv(ftests_Mcali, "output/tables/ftests_Mcali.csv", row.names=F)

# Prep for graphing
ftests_Mcali$class_graphing <- ftests_Mcali$class
ftests_Mcali$class_graphing <- str_remove_all(ftests_Mcali$class_graphing, "_variant")
ftests_Mcali$class_graphing <- gsub("_", " ", ftests_Mcali$class_graphing)

# Order classes based on OR
ftests_Mcali <- ftests_Mcali %>% mutate(class_graphing = fct_reorder(class_graphing, OR))

# Graph OR of Fishers exact tests
pdf("output/figures/outlier_analyses/Fishers_exact_test_Mcali.pdf", width = 8, height = 8)
ggplot(ftests_Mcali, aes(x = log2(OR), y = class_graphing, fill = -log(p.fet))) + 
  geom_vline(xintercept = 0, col="black", linetype="dashed") + ylab("") + 
  geom_linerange(aes(xmin = log2(lci), xmax = log2(uci)), linewidth = 1) +
  geom_point(shape = 21, size = 8) +
  scale_fill_gradient(low = "#e6e4e4", high = "#b55c04", name="-log10(p)") + 
  #labs(title = expression(paste(italic("M. californianus"), " cross-sectional thickness"))) +
  theme_linedraw(base_size=30) + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.margin = margin(t = 10, r = 50, b = 10, l = 10,, unit = "pt")) +
  theme(legend.title = element_text(size = 19), legend.text = element_text(size = 20), legend.position = c(0.218, 0.83), legend.background = element_rect(color = "black", fill = "white", linewidth = 0.5, linetype = "solid"))
dev.off()

