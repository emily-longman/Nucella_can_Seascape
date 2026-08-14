# Graph Baypass XtX output.

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
source("/gpfs1/home/e/l/elongman/software/baypass_public/utils/baypass_utils.R")
install.packages(c('data.table', 'dplyr', 'ggplot2', 'mvtnorm', 'geigen'))
library(data.table)
library(dplyr)
library(ggplot2)
library(mvtnorm)
library(geigen)
library(tidyverse)
library(foreach)

# ================================================================================== #

# Read in Baypass results
NC.omega <- as.matrix(read.table("data/processed/GEA/baypass/omega/NC_baypass_mat_omega.out"))
XtX <- read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_pi_xtx.out", header=T)
pops <- read.table("data/processed/fastq_to_vcf/guide_files/N.canaliculata_pops.vcf_pop_names.txt", header=F)
snp.meta <- read.table("data/processed/GEA/baypass/snpdet", header=F)

# Re-name pool names 
colnames(NC.omega) <- pops$V1
rownames(NC.omega) <- pops$V1

# Re-name snp metadata
colnames(snp.meta) <- c("chr", "pos", "allele1", "allele2")
 
# ================================================================================== #

# Check the behavior of the p-values associated to the XtXst estimator
pdf("output/figures/GEA/Baypass_xtx_hist.pdf", width = 5, height = 5)
hist(10**(-1*XtX$log10.1.pval.), freq=F, breaks=50)
abline(h=1, col="red")
dev.off()

# Graph xtx
pdf("output/figures/GEA/Baypass_xtx.pdf", width = 5, height = 5)
plot(XtX$XtXst)
dev.off()

# Graph outliers
pdf("output/figures/GEA/Baypass_xtx_outliers.pdf", width = 5, height = 5)
plot(XtX$log10.1.pval., ylab="XtX P-value (-log10 scale)" )
abline(h=3, lty=2, col="red") #0.001 p-value threshold
dev.off()



# ================================================================================== #


# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# Rank normalization

L = dim(SNP.XtX.dt)[1] 
    
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

# Then do an alpha of 0.01 or 0.05

# ================================================================================== #
# ================================================================================== #

# Compare POD and original data estimates

# Get estimate of omega from POD analysis
NC.POD.omega <- as.matrix(read.table("data/processed/GEA/baypass/POD/NC_baypass_POD_mat_omega.out"))

# Compare POD omega and original omega
pdf("output/figures/GEA/Baypass_compare_omega_POD.pdf", width = 5, height = 5)
plot(NC.POD.omega, NC.omega)
abline(a=0, b=1)
dev.off()
# Note: you want similar values between the two .... not the case with the Nucella data

# Get the Forstner and Moonen Distance (FMD) between simulated and original posterior estimates (here a smaller value is better) 
fmd.dist(NC.POD.omega, NC.omega) # 5.374831 --- this seems quite high...

# Get estimates (posterior mean) of both the a_pi and b_pi parameters of the Pi Beta distribution from the POD analysis
pi.beta.coef=read.table("data/processed/GEA/baypass/xtx/NC_baypass_core_summary_beta_params.out", h=T)$Mean
pod.pi.beta.coef=read.table("data/processed/GEA/baypass/POD/NC_baypass_POD_summary_beta_params.out", h=T)$Mean

# Graph the estimate
pdf("output/figures/GEA/Baypass_compare_pi_beta_coef_POD.pdf", width = 5, height = 5)
plot(pod.pi.beta.coef, pi.beta.coef)
abline(a=0, b=1)
dev.off()

# ================================================================================== #

# XtX calibration

# Get the POD XtX
POD.XtX=read.table("data/processed/GEA/baypass/POD/NC_baypass_POD_summary_pi_xtx.out")$M_XtX

# Compute the 1% threshold (i.e., identify SNPs where the xtx values are above the 99% significance threshold from the POD)
pod.thres=quantile(POD.XtX, probs=0.99) #NA

# Add the threshold to the actual XtX plot
pdf("output/figures/GEA/Baypass_xtx_POD_thres.pdf", width = 5, height = 5)
plot(XtX$XtXst)
abline(h=pod.thres, lty=2) 
dev.off()

# ... Problem.... getting NA

# ================================================================================== #
# ================================================================================== #

# Merge baypass results and SNP metadata 
SNP.XtX <- cbind(snp.meta, XtX)
SNP.XtX.dt <- as.data.table(SNP.XtX)

# ================================================================================== #

# Import local score functions (https://forge-dga.jouy.inra.fr/documents/809).
source('data/processed/GEA/baypass/xtx/localscore/scorelocalfunctions.R')

# Create key based on the chromosome
setkey(SNP.XtX.dt, chr)
# Total number of unique chromosomes: 18,369
Nchr=length(SNP.XtX.dt[,unique(chr)])

# Compute the absolute position in the genome. This is useful for doing genome-wide plots. 
# Add this column (posT) to the data table. 
chrInfo=SNP.XtX.dt[,.(L=.N, cor=autocor(log10.1.pval.)), chr]
setkey(chrInfo, chr)
tmp=data.table(chr=SNP.XtX.dt[,unique(chr),], S=cumsum(c(0,chrInfo$L[-Nchr])))
setkey(tmp,chr)
SNP.XtX.dt[tmp,posT:=pos+S]


# Choice of $\xi$ (1,2,3 or 4)

# To choose the apropiate threshold ($\xi = 1$ or 2) we look at the distribution of −log10(p − value). Then the score function will be $−log10(p − value) − \xi$. 

pdf("output/figures/GEA/Baypass_xtx_lindley_qplot.pdf", width = 5, height = 5)
qplot(log10.1.pval., data=SNP.XtX.dt, geom='histogram', binwidth=0.1, main='P-values histogram')
dev.off()

# Remember that we have to choose some value between mean(-log10(mydata$pval)) and max(-log10(mydata$pval)).
mean(SNP.XtX.dt$log10.1.pval.)  #0.4177791
max(SNP.XtX.dt$log10.1.pval.) #10.44488

#pdf("output/figures/GEA/Baypass_xtx_lindley_plot.pdf", width = 5, height = 5)
#p<-ggplot(data=SNP.XtX.dt, aes(posT,log10.1.pval.)) + geom_point()
#p + geom_abline(intercept=c(1,2), slope=0) 
#dev.off()

# Computation of the score and the Lindley Process

xi=1
SNP.XtX.dt[,score:= log10.1.pval.-xi]
# The score mean must be negative
mean(SNP.XtX.dt$score) # -1.582221
SNP.XtX.dt[,lindley:=lindley(score),chr]

# Compute significance threshold for each chromosome

# Uniform distribution of p-values

# If the distribution of the p-values is uniform and if $\xi=1$ or 2, 
# it is possible to compute the thresholds of significancy for each chromosome directely, given the length and the autocorrelation of the chromosome .

chrInfo[,th:=thresUnif(L, cor, 2),chr]
(SNP.XtX.dt=SNP.XtX.dt[chrInfo])
(sigZones=SNP.XtX.dt[chrInfo, sig_sl(lindley, pos, unique(th)),chr])


## If the distribution of the p-values is not uniform, then a re-sampling strategy should be used. 
# Note: needed to edit 'scorelocalfunctions.R' to switch 'p-val' to 'log10.1.pval.'

setwd("output/figures/GEA")

coefsG=coefsGumb(SNP.XtX.dt, Ls=seq(10000,70000,10000), nSeq=5000) # this may take a while.
# Many of the coefficients are NA .....so can't calculate a threshold

chrInfo[,thG05:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.05),]
chrInfo[,thG01:=threshold(L, cor, coefsG$aCoef, coefsG$bCoef,0.01),]

SNP.XtX.dt=SNP.XtX.dt[chrInfo]

sigZones05=SNP.XtX.dt[,sig_sl(lindley, pos, unique(thG05)),chr]
ind=which(sigZones05[,peak]>0)
write.table(sigZones05[ind,],file='SL_xi1_signif5.txt',col.names=T,row.names=F,quote=F)

sigZones01=SNP.XtX.dt[,sig_sl(lindley, pos, unique(thG01)),chr]
ind=which(sigZones01[,peak]>0)
write.table(sigZones01[ind,],file='SL_xi1_signif1.txt',col.names=T,row.names=F,quote=F)
