# Nucella ABC analysis - ABC estimation

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'foreach', 'lme4', 'abc', 'RColorBrewer', 'minpack.lm'))
library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(foreach)
library(lme4)
library(abc)
library(RColorBrewer)
library(minpack.lm)

# ================================================================================== #

# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))[-4]

# ================================================================================== #

# Generate output directories

# Data directory
out_data_dir <- paste("data/processed/SLiM/Mcali_ABC")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/Mcali_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T)
# Extract just sites and env var
Mcali <- afs.Mcali[,c(7,13)] %>% distinct()

# Load data
real_All <- get(load("data/processed/SLiM/Mcali_ABC/real_data.Rdata"))
#sim_All <- get(load("data/processed/SLiM/Mcali_ABC/sim_data.Rdata"))
#sim_All <- get(load("data/processed/SLiM/Mcali_ABC/sim_data_v2.Rdata"))
sim_All <- get(load("data/processed/SLiM/Mcali_ABC/sim_data_v3.Rdata"))
sim_All <- get(load("data/processed/SLiM/Mcali_ABC/sim_data_v4.Rdata"))

# ================================================================================== #

# Estimate means
real_data = dplyr::select(ungroup(real_All), 
                          #Fsg, Fgt, Fst,
                          Fst.thin.mod, Fst.mod.thick, Fst.thin.thick,
                          deltaAF.thick.mod, deltaAF.thick.thin,
                          cor.pearson, 
                          cor.spearman, 
                          corAF.pearson, 
                          corAF.spearman, 
                          #fix1, 
                          #fix0, 
                          #poly, 
                          mean.AF,
                          #asym, xmid, scal,
                          p0, p1, p2, p4, p5, p6,    
                          p7, p8, p9, p10, p11, p12,
                          p13, p14, p15, p16, p17, p18
                          )
simulated_data = dplyr::select(ungroup(sim_All), 
                          #Fsg, Fgt, Fst, 
                          Fst.thin.mod, Fst.mod.thick, Fst.thin.thick,
                          deltaAF.thick.mod, deltaAF.thick.thin,
                          cor.pearson, 
                          cor.spearman, 
                          corAF.pearson, 
                          corAF.spearman, 
                          #fix1, 
                          #fix0, 
                          #poly, 
                          mean.AF, 
                          p0, p1, p2, p4, p5, p6,    
                          p7, p8, p9, p10, p11, p12,
                          p13, p14, p15, p16, p17, p18
                          )
sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), thresh, k, m), as.numeric)) 
#sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), m, thresh, k, mag, N ), as.numeric)) 

# ================================================================================== #

# Do ABC
# When method is "loclinear", a local linear regression method corrects for the imperfect match between S(y) and S(y0).
abc_fit <- abc(
  target = real_data,
  param = sim_parameters,
  sumstat = simulated_data,
  tol = 0.1,
  method = "loclinear")

# ================================================================================== #

# Extract posterior
post <- as.data.frame(abc_fit$adj.values)
# Reformat
post_long <- pivot_longer(
  post,
  cols = everything(),
  names_to = "parameter",
  values_to = "value")

# Graph posteriors
pdf("output/figures/SLiM/Mcali_ABC/Mcali_posteriors_v4_fewerstats.pdf", width = 5, height = 8)
ggplot(post_long, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(parameter~. , scales = "free", ncol = 1) +
  theme_bw(base_size=20)
dev.off()

# ================================================================================== #

# Report the simulation that best match the data

# Mean posterior for each parameter
post_long %>%
  group_by(parameter) %>%
  summarise(mean.par = mean(value, na.rm = T))

# Extract best parameters
best <- which.max(abc_fit$weights)
abc_fit$unadj.values[best, ]

# Create var with best value for each param
m_best = abc_fit$unadj.values[best,"m"]
k_best = abc_fit$unadj.values[best,"k"]
thresh_best = abc_fit$unadj.values[best,"thresh"]
#mag_best = abc_fit$unadj.values[best,"mag"]
#N_best = abc_fit$unadj.values[best,"N"]

# Assigned Mcali values
env = c(
1.801197107,
1.758352318,
1.715225361,
1.677236687,
1.734608112,
1.844818925,
1.735326540,
1.615483995,
1.961925234,
1.738814502,
1.517563421,
1.487237983,
1.567637305,
1.472519437,
1.595566753,
2.206475463,
1.520675474,
1.431920910);


# ================================================================================== #
# ================================================================================== #

# Compare to real data

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Allele frequency data for top ph hits
afs.Mcali <- read.csv("data/processed/SLiM/afs.McaliThk.outlier.csv", header=T) %>% mutate(nsnails = 20)
# Order by site/latitude
afs.Mcali %<>% arrange(desc(latitude))
# Merge with ecovar
afs.Mcali.sims.tmp <- cbind(afs.Mcali, env)
# Join
afs.Mcali.sims <- left_join(afs.Mcali.sims.tmp, dplyr::select(ecovars, "Site", "sim_eq", "Demographic Cluster"))


# Load top SNP
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595") %>%
  left_join(dplyr::select(ecovars, Site, ph_mean))

# Reformat Sim data
# Extract p0-p18
#sim_sub <- sim_All[,c(12:35)]
sim_sub <- sim_All[,c(17:40)]
# Change to long format
#sim_Data.melt <- sim_sub %>%
#  reshape2::melt(id = c("repId","m","thresh","k_1","k_2","N"),
#                 variable.name = "sim_eq",
#                 value.name = "AF_true")

sim_Data.melt <- sim_sub %>%
  reshape2::melt(id = c("repId","m","thresh","k","mag","N"),
                 variable.name = "sim_eq",
                 value.name = "AF_true")

#sim_All[which(sim_All$m == m_best & sim_All$thresh == thresh_best & sim_All$k_1 == k_1_best & sim_All$k_2 == k_2_best),]
#sim_All[which(sim_All$m == m_best & sim_All$thresh == thresh_best & sim_All$k == k_best & sim_All$mag == mag_best & sim_All$N == N_best),]
sim_All[which(sim_All$thresh == thresh_best & sim_All$k == k_best),]

# ================================================================================== #

# Graph real data 

# Make Site factor
afs.Mcali.sims$Site <- factor(afs.Mcali.sims$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Graph real
pdf("output/figures/SLiM/Mcali_ABC/real_AFs_topsnp.pdf", width = 5, height = 5)
ggplot(afs.Mcali.sims, aes(x = AF, y = mean_integrated_thk, fill = Site)) + geom_point(size = 5, shape = 21) +  scale_fill_manual(values = mycolors) + 
  labs(x= "AF", y= "Shell Thk") +
  theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()

####

# Graph best sim

# Extract best
sim_AFs_best <- sim_Data.melt %>%
  filter(thresh == thresh_best, k == k_best) %>%
  left_join(dplyr::select(ecovars, Site, sim_eq, shell_thk))
# Make Site factor
sim_AFs_best$Site <- factor(sim_AFs_best$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))


# Graph sim
pdf("output/figures/SLiM/Mcali_ABC/sim_AFs_best_v4.pdf", width = 5, height = 5)
ggplot(sim_AFs_best, aes(x = AF_true, y = shell_thk, fill = Site)) + geom_point(size = 5, shape = 21) + scale_fill_manual(values = mycolors) +
labs(x= "AF", y= "Shell Thk") +
theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()

# Calculate the AF average of the iterations of the best param
sim_AFs_best_mean <- sim_AFs_best %>%
  group_by(m, thresh, k, mag, N, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, shell_thk))
  # Make Site factor
sim_AFs_best_mean$Site <- factor(sim_AFs_best_mean$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/Mcali_ABC/sim_AFs_best_repmeans_v4.pdf", width = 5, height = 5)
ggplot(sim_AFs_best_mean, aes(x = AF_true_mean, y = shell_thk, fill = Site)) + geom_point(size = 5, shape = 21) + scale_fill_manual(values = mycolors) + 
labs(x= "AF", y= "Shell Thk") +
theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()


# Fit using self-starting parameters - getting error
#sim_AFs_best_mean_sub <- sim_AFs_best_mean[,c("AF_true_mean","shell_thk")]
#mod_sim <- nlsLM(AF_true_mean ~ SSlogis(shell_thk, Asym, xmid, scal), data = sim_AFs_best_mean_sub)
#mod_sim_fit <- coef(mod_sim)


# Selection
sel <- function(x, z, k) {
  1 / (1 + exp((x - z)/k)) - 0.5
}

# Graph
pdf("output/figures/SLiM/Mcali_ABC/Mcali_selection.pdf", width = 5, height = 5)
plot(afs.Mcali$mean_integrated_thk, sel(afs.Mcali$mean_integrated_thk, 1.8, 4))
dev.off()
pdf("output/figures/SLiM/Mcali_ABC/Mcali_selection_zoomed_out.pdf", width = 5, height = 5)
plot(seq(0, 5, by = 0.01), sel(seq(0, 5, by = 0.01), 1.79, 4))
dev.off()


##############

# Visually look at sims to see best match
# Extract best
sim_AFs_test <- sim_Data.melt %>%
  #filter(thresh == "1.77", k == "0.1", m == "0.0001") %>%
  filter(thresh == "2", k == "0.2") %>%
  left_join(dplyr::select(ecovars, Site, sim_eq, shell_thk))
# Make Site factor
sim_AFs_test$Site <- factor(sim_AFs_test$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/Mcali_ABC/sim_AFs_best_v1_test.pdf", width = 5, height = 5)
ggplot(sim_AFs_test, aes(x = AF_true, y = shell_thk, fill = Site)) + geom_point(size = 5, shape = 21) + scale_fill_manual(values = mycolors) +
labs(x= "AF", y= "Shell Thk") +
theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()

# Calculate the AF average of the iterations of the best param
sim_AFs_test_mean <- sim_AFs_test %>%
  group_by(m, thresh, k, mag, N, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, shell_thk))
  # Make Site factor
sim_AFs_test_mean$Site <- factor(sim_AFs_test_mean$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/Mcali_ABC/sim_AFs_best_repmeans_v1_test.pdf", width = 5, height = 5)
ggplot(sim_AFs_test_mean, aes(x = AF_true_mean, y = shell_thk, fill = Site)) + geom_point(size = 5, shape = 21) + scale_fill_manual(values = mycolors) + 
labs(x= "AF", y= "Shell Thk") +
theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()
