# Analyze omega matrix

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
install.packages(c('corrplot', 'ggplot2', 'RColorBrewer', 'WriteXLS', 'ape'))
library(corrplot)
library(ggplot2)
library(RColorBrewer)
library(WriteXLS)
library(ape)

# ================================================================================== #

# Load data
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
NC.omega <- as.matrix(read.table("data/processed/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out"))
XtX <- read.table("data/processed/outlier_analyses/baypass/omega/NC_baypass_summary_pi_xtx.out", header=T)

# ================================================================================== #

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

write.table(NC.omega, "output/tables/NC.omega.txt", sep='\t')

# Re-order populations so in latitudinal order
NC.omega_reorder <- as.matrix(read.table("output/tables/NC.omega_reorder.txt", header=T))

# ================================================================================== #

# Create a correlation matrix of the omega values -- assess genomic differentiation between pools
cor.mat <- cov2cor(NC.omega_reorder)

# Graph correlation matrix
pdf("output/figures/outlier_analyses/Omega/Baypass_omega_cor_matrix_title.pdf", width = 5, height = 5)
corrplot(cor.mat, method = "color", mar=c(2,1,2,2)+0.1, main=expression("Correlation map based on"~hat(Omega)))
dev.off()
# Graph correlation matrix (no title)
pdf("output/figures/outlier_analyses/Omega/Baypass_omega_cor_matrix.pdf", width = 5, height = 5)
corrplot(cor.mat, method = "color", mar=c(2,1,2,2)+0.1)
dev.off()

# ================================================================================== #

# Assess population differentiation with hierarchical clustering
NC.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))

# Graph tree 
pdf("output/figures/outlier_analyses/Omega/Baypass_hier_clustering_title.pdf", width = 5, height = 5)
plot(NC.tree,type="p",
     main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))
dev.off()

# Graph tree (no title)
pdf("output/figures/outlier_analyses/Omega/Baypass_hier_clustering.pdf", width = 5, height = 5)
plot(NC.tree,type="p")
dev.off()

# ================================================================================== #
# ================================================================================== #

# Graph estimates of the XtX differentation measure

# Check the behavior of the p-values associated to the XtXst estimator
pdf("output/figures/outlier_analyses/Omega/Baypass_xtx_hist_omega.pdf", width = 5, height = 5)
hist(10**(-1*XtX$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/outlier_analyses/Omega/Baypass_xtx_omega.pdf", width = 5, height = 5)
plot(XtX$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/outlier_analyses/Omega/Baypass_xtx_outliers_omega.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="-log10(XtX P-value)")
abline(h=3, lty=2, col="red") #0.001 p-value threshold
dev.off()