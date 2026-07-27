# Nucella ABC analysis - ABC estimation

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
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# ================================================================================== #

# Generate output directories

# Data directory
out_data_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

# Figure directory
out_fig_dir <- paste("output/figures/SLiM/ph_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data
afs.ph <- read.csv("data/processed/SLiM/afs.ph.g27343.BF.POD.csv", header=T)
ph <- afs.ph[, c(2,8)] %>% distinct()

# Load data
real_All <- get(load("data/processed/SLiM/ph_ABC/real_data.Rdata"))
#sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars2.Rdata"))
#sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars3.Rdata"))
#sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars4.Rdata"))
#sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars5.Rdata"))

sim_All <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars6.Rdata"))
#sim_6 <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars6.Rdata"))
#sim_6_mSet <- get(load("data/processed/SLiM/ph_ABC/sim_data_morevars6_mSet.Rdata"))

#sim_All <- rbind(sim_6[which(sim_6$m=="0.001"),], sim_6_mSet)
# ================================================================================== #

# Estimate means

real_data = dplyr::select(ungroup(real_All), 
                          #Fsg, Fgt, Fst, 
                          Fst.L.M, Fst.M.H, Fst.L.H,
                          deltaAF.L.M, deltaAF.H.M,
                          cor.pearson, cor.spearman, 
                          corAF.pearson, corAF.spearman, 
                          #fix1, fix0, poly,
                          mean.AF,
                          #asym, xmid, scal,
                          #p0, p1, p2, p3, p4, 
                          #p5, #p6,
                          #p7, p8, p9, p10, p11, p12,
                          #p13, p14, p15, p16, p17, p18
                          )
simulated_data = dplyr::select(ungroup(sim_All), 
                          #Fsg, Fgt, Fst, 
                          Fst.L.M, Fst.M.H, Fst.L.H,
                          deltaAF.L.M, deltaAF.H.M, 
                          cor.pearson, cor.spearman, 
                          corAF.pearson, corAF.spearman, 
                          #fix1, fix0, poly, 
                          mean.AF,
                          #p0, p1, p2, p3, p4, 
                          #p5, #p6,
                          #p7, p8, p9, p10, p11, p12,
                          #p13, p14, p15, p16, p17, p18
                          )

#sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), thresh, k), as.numeric)) 
sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), thresh, k, m), as.numeric)) 
#sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), m, thresh, k, mag, N ), as.numeric)) 
#sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), thresh, k_1, k_2 ), as.numeric)) # EXCLUDE m and N as those are invariant

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
pdf("output/figures/SLiM/ph_ABC/pH_posteriors_morvars6.pdf", width = 5, height = 5)
ggplot(post_long, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(parameter~. , scales = "free", ncol = 1) +
  theme_bw(base_size=15)
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
#kb2 = abc_fit$unadj.values[best,"k_2"]
thresh_best = abc_fit$unadj.values[best,"thresh"]
#mag_best = abc_fit$unadj.values[best,"mag"]
#N_best = abc_fit$unadj.values[best,"N"]

### plots
env = c(
  8.023377561,
  8.016675575,
  8.013330447,
  8.006025013,
  7.987332623,
  7.926713508,
  7.957538759,
  7.940880291,
  7.937233951,
  7.933377415,
  7.923570354,
  7.941466321,
  7.941527015,
  7.946801475,
  7.963869108,
  7.96609182,
  7.948259634,
  7.946165852,
  7.950069857    
);



sel <- function(x, z, k) {
  1 / (1 + exp((x - z)/k)) - 0.5
}
# Graph
pdf("output/figures/SLiM/ph_ABC/pH_selection_curve.pdf", width = 5, height = 5)
plot(ph$ph_mean, sel(ph$ph_mean, 7.985, 0.04))
dev.off()
pdf("output/figures/SLiM/ph_ABC/pH_selection_curve_biggerph.pdf", width = 5, height = 5)
plot(seq(7.8, 8.1, by = 0.01), sel(seq(7.8, 8.1, by = 0.01), 7.988, 0.12))
dev.off()

# ================================================================================== #
# ================================================================================== #

# Compare to real data

# Load ecovars
ecovars <- fread("guide_files/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

# Load top SNP
phafs <- fread("data/processed/baypass/afs.ph.g27343.BF.POD.csv") %>% mutate(nsnails = 20)
# Extract top pH hit
topsnp <- phafs %>%
  filter(SNP_id == "ntLink_3821_1595") %>%
  left_join(dplyr::select(ecovars, Site, ph_mean))

# Reformat Sim data
# Extract p0-p18
sim_sub <- sim_All[,c(17:41)]
#sim_sub <- sim_All[,c(12:36)]
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
#sim_All[which(sim_All$thresh == thresh_best & sim_All$k == k_best & sim_All$m == m_best),]
sim_All[which(sim_All$thresh == thresh_best & sim_All$k == k_best),]

# ================================================================================== #

# Graph Real

# Make Site factor
topsnp$Site <- factor(topsnp$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
# Graph real
pdf("output/figures/SLiM/ph_ABC/real_AFs_topsnp.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = AF, y = ph_mean, fill = Site)) + geom_point(size = 6, shape = 21) +  
  scale_fill_manual(values = mycolors) + labs(x="AF", y="Mean pH") + theme_linedraw(base_size=24) + theme(legend.position = "none")
dev.off()

####

# Graph Best Sim

# Extract best
sim_AFs_best <- sim_Data.melt %>%
  #filter(thresh == thresh_best, k == k_best) %>%
  filter(thresh == thresh_best, k == k_best, m == m_best) %>%
  #filter(m == m_best, thresh == thresh_best, k == k_best, mag == mag_best, N == N_best) %>%
  #filter(m == mig_best, thresh == thresh_best, k_1 == k_1_best, k_2 == k_2_best) %>% 
  #filter(thresh == thresh_best, k_1 == k_1_best, k_2 == k_2_best) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, ph_mean))
# Make Site factor
sim_AFs_best$Site <- factor(sim_AFs_best$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/ph_ABC/sim_AFs_best_morevars6.pdf", width = 5, height = 5)
ggplot(sim_AFs_best, aes(x = AF_true, y = ph_mean, fill = Site)) + geom_point(size = 5, shape = 21, alpha=0.8) + 
  scale_fill_manual(values = mycolors) + labs(x="AF", y="Mean pH") + theme_linedraw(base_size=24) + theme(legend.position = "none")
dev.off()

# Calculate the AF average of the iterations of the best param
sim_AFs_best_mean <- sim_AFs_best %>%
  #group_by(m, thresh, k_1, k_2, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  group_by(m, thresh, k, mag, N, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, ph_mean))
  # Make Site factor
sim_AFs_best_mean$Site <- factor(sim_AFs_best_mean$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/ph_ABC/sim_AFs_best_repmeans_morevars6.pdf", width = 5, height = 5)
ggplot(sim_AFs_best_mean, aes(x = AF_true_mean, y = ph_mean, fill = Site)) + geom_point(size = 5, shape = 21) + 
  scale_fill_manual(values = mycolors) + labs(x="AF", y="Mean pH") + theme_linedraw(base_size=24) + theme(legend.position = "none")
dev.off()


# Fit sigmoid
sim_AFs_best_mean_sub <- sim_AFs_best_mean[,c("AF_true_mean","ph_mean")]
mod_sim <- nlsLM(AF_true_mean ~ SSlogis(ph_mean, Asym, xmid, scal), data = sim_AFs_best_mean_sub)
mod_sim_fit <- coef(mod_sim)

##############

# Visually look at sims to see best match
# Extract best
sim_AFs_test <- sim_Data.melt %>%
  filter(thresh == 7.99, k == 0.07, m == 0.001) %>%
  left_join(dplyr::select(ecovars, Site, sim_eq, ph_mean))
# Make Site factor
sim_AFs_test$Site <- factor(sim_AFs_test$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))

# Graph sim
pdf("output/figures/SLiM/ph_ABC/sim_AFs_best_v6_test.pdf", width = 5, height = 5)
ggplot(sim_AFs_test, aes(x = AF_true, y = ph_mean, fill = Site)) + geom_point(size = 5, shape = 21) + 
  scale_fill_manual(values = mycolors) + labs(x= "AF", y= "Mean pH") + theme_linedraw(base_size=26) + theme(legend.position = "none")
dev.off()

sim_AFs_test_mean <- sim_AFs_test %>%
  #group_by(m, thresh, k_1, k_2, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  group_by(m, thresh, k, mag, N, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, ph_mean))
# Fit sigmoid
sim_AFs_test_mean_sub <- sim_AFs_test_mean[,c("AF_true_mean","ph_mean")]
mod_sim_test <- nlsLM(AF_true_mean ~ SSlogis(ph_mean, Asym, xmid, scal), data = sim_AFs_test_mean_sub)
mod_sim_test_fit <- coef(mod_sim_test)
