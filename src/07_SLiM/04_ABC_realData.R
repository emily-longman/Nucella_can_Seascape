# Nucella ABC analysis
## Part 1 - real data

# Clear memory
rm(list=ls())

# Stop exponential
options(scipen = 999)

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'lme4'))
library(data.table)
library(tidyverse)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(lme4)

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load and format data

# Ecological variables
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Allele frequency data for top ph hits
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit and join with eco vars
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595")  %>%
  left_join(dplyr::select(ecovars, Site, Latitude, sim_eq)) %>%
  arrange(-Latitude)

# Make a pooled object with the raw data
pool.real <- new("pooldata",
                 npools=19, #### Rows = Number of pools
                 nsnp=1, ### Columns = Number of SNPs
                 refallele.readcount=as.matrix(t(topsnp$Count)),
                 readcoverage=as.matrix(t(topsnp$Cov)),
                 poolsizes=topsnp$nsnails * 2,
                 poolnames = topsnp$Site)

# ================================================================================== #

# Estimate the key statistics

# Hierach fst between N vs S
fst.phylogeo <- computeFST(pool.real,
                        method = "Anova",
                        struct = topsnp$shape, verbose = FALSE)
# "FST": estimate of genome-wide Fst over all the populations
# "FSG": estimate of genome-wide within-group differentiation (Fsg)
# "FGT": estimate of genome-wide between-group differentiation (Fgt)

# Raw correlation between mean pH and AF
rawcor = cor.test(~ ph_mean+AF, data = topsnp)

# Correlation between AF and AF - 1 for real data
rawcorAF = cor.test(topsnp$AF, topsnp$AF)

# ================================================================================== #

# Extract AFs of top SNP and format similar to sim output
AF_pool = data.frame(t(topsnp$AF))
names(AF_pool) = topsnp$sim_eq

# Format Data
real_data =
data.frame(
  Fsg = fst.phylogeo$snp.Fstats[1],
  Fgt = fst.phylogeo$snp.Fstats[2],
  Fst = fst.phylogeo$snp.Fstats[3],
  cor = rawcor$estimate,
  corAF = rawcorAF$estimate,
  fix1 = sum(topsnp$AF ==1),
  fix0 = sum(topsnp$AF ==0),
  poly = sum(topsnp$AF  > 0 & topsnp$AF  < 1),
  mean.AF = mean(topsnp$AF),
  AF_pool)

# ================================================================================== #

# Save output
save(real_data, file = "data/processed/SLiM/ph_ABC/real_data.Rdata")





# ================================================================================== #
# ================================================================================== #

library(stats)

# Make Site factor
topsnp$Site <- factor(topsnp$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Graph real
pdf("output/figures/SLiM/real_AFs_topsnp_sigmoid.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = AF, y = ph_mean, fill = Site)) + geom_point(size = 3, shape = 21) +  scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()
pdf("output/figures/SLiM/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = ph_mean, y = AF, fill = Site)) + geom_point(size = 3, shape = 21) +  scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()

# Subset data
topsnp_sub <- topsnp[,c(5,8)]

# Fit using self-starting parameters
mod <- nls(AF ~ SSlogis(ph_mean, Asym, xmid, scal), data = topsnp_sub)
mod_fit <- coef(mod)

pdf("output/figures/SLiM/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
plot(topsnp_sub$ph_mean, topsnp_sub$AF, pch = 20)
curve(SSlogis(x, mod_fit["Asym"], mod_fit["xmid"], mod_fit["scal"]), lwd = 2, col = 'lightblue', add = TRUE)
dev.off()

pdf("output/figures/SLiM/real_AFs_topsnp_sigmoid_flipped.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = ph_mean, y = AF, fill = Site)) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = mycolors) +
stat_function(fun = SSlogis, args = list(Asym = mod_fit["Asym"], xmid = mod_fit["xmid"], scal = mod_fit["scal"])) + theme_linedraw()
dev.off()



# Other ideas
# Selection
s <- function(x, asym, thresh, scal) {
  -1 * (asym/(1 + exp((x-thresh)/scal)) - asym/2)
}
# Load data
afs.ph <- read.csv("data/processed/SLiM/afs.ph.g27343.BF.POD.csv", header=T)
ph <- afs.ph[, c(2,8)] %>% distinct()

# Graph
pdf("output/figures/SLiM/pH_dist_sel_testing.pdf", width = 5, height = 5)
plot(ph$ph_mean, s(ph$ph_mean, 0.986, 7.98, -0.0043))
dev.off()
