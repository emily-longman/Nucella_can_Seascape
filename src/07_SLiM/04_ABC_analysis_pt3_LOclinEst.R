# Nucella ABC analysis
## Part 3 - ABC estimation

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
#install.packages(c('data.table', 'tidyverse', 'magrittr', 'reshape2', 'gmodels', 'poolfstat', 'foreach', 'lme4', 'abc', 'RColorBrewer'))
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

# ================================================================================== #

# Generate output directories

# Figure directory
out_fig_dir <- paste("data/processed/SLiM/ph_ABC")
if (!dir.exists(out_fig_dir)) {dir.create(out_fig_dir)}

# ================================================================================== #

# Load data

real_All <- get(load("data/processed/SLiM/ph_ABC/real.data.Rdata"))
sim_All <- get(load("data/processed/SLiM/ph_ABC/sim.results.Rdata"))

# ================================================================================== #

# Estimate means

real_data = dplyr::select(ungroup(real_All), 
                          Fsg, Fgt, Fst, 
                          cor, #fix1, fix0, poly, 
                          #mean.AF,
                          #p0, p1, p2, p3, p4, p5, p6,    
                          #p7, p8, p9, p10, p11, p12,
                          #p13, p14, p15, p16, p17, p18
                          )
simulated_data = dplyr::select(ungroup(sim_All), 
                               Fsg, Fgt, Fst, 
                               cor, #fix1, fix0, poly, 
                               #mean.AF,
                               #p0, p1, p2, p3, p4, p5, p6,    
                               #p7, p8, p9, p10, p11, p12,
                               #p13, p14, p15, p16, p17, p18
                               )
sim_parameters = as.data.frame(lapply(dplyr::select(ungroup(sim_All), thresh,   k_1,   k_2 ), as.numeric))
# EXCLUDE m and N as those are invariant

# ================================================================================== #

# Do ABC
# When method is "loclinear", a local linear regression method corrects for the imperfect match between S(y) and S(y0).
abc_fit <- abc(
  target = real_data,
  param = sim_parameters,
  sumstat = simulated_data,
  tol = 0.1,
  method = "loclinear"
)

# ================================================================================== #

# Extract posterior
post <- as.data.frame(abc_fit$adj.values)
post_long <- pivot_longer(
  post,
  cols = everything(),
  names_to = "parameter",
  values_to = "value")

# Graph posteriors
pdf("output/figures/SLiM/pH_posteriors.pdf", width = 5, height = 5)
ggplot(post_long, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(parameter~. , scales = "free", ncol = 1) +
  theme_bw()
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
kb1 = abc_fit$unadj.values[best,"k_1"]
kb2 = abc_fit$unadj.values[best,"k_2"]
thre = abc_fit$unadj.values[best,"thresh"]

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


# Logistic sigmoid scaled to [-1, 1]
f <- function(x, z, k) {
  2/(1+exp(-(x - z)/k))-1
}

# Plot sigmoid of best fit param
pdf("output/figures/SLiM/sel_curve.pdf", width = 5, height = 5)
data.frame(env=env, s=f(env, thre, kb1)) %>%
  ggplot(aes(x=env, y=s)) + geom_line(linewidth=2) + theme_linedraw()
dev.off()

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


####
# Load sim data
sim_Data <- get(load("data/processed/SLiM/ph_results.Rdata"))
names(sim_Data)[1:19] = paste("p", 0:18, sep ="")
# Reformat
sim_Data.melt <- sim_Data %>%
  reshape2::melt(id = c("repId","m","thresh","k_1","k_2",
                        "N","state","sim.cycle"),
                 variable.name = "sim_eq",
                 value.name = "AF_true")

# ================================================================================== #

###recall...
kb1 = abc_fit$unadj.values[best,"k_1"]
kb2 = abc_fit$unadj.values[best,"k_2"]
thre = abc_fit$unadj.values[best,"thresh"]

# What is w? 
w <- abc_fit$weights / sum(abc_fit$weights)
colSums(abc_fit$unadj.values * w)


# Plot
pdf("output/figures/SLiM/sim_AFs.pdf", width = 5, height = 5)
sim_Data.melt %>% 
  filter( 
          thresh == thre,
          k_1 == kb1,
          k_2 ==kb2 
          ) %>%
  left_join(dplyr::select(ecovars, sim_eq, ph_mean)) %>%
  ggplot(aes(y=ph_mean, x=1-AF_true)) + geom_point(size = 2.0) 
dev.off()



#####
# Color palette
nb.cols <- 19
mycolors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(nb.cols))

# Graph real
topsnp$Site <- factor(topsnp$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
pdf("output/figures/SLiM/real_AFs_topsnp.pdf", width = 5, height = 5)
ggplot(topsnp, aes(x = AF, y = ph_mean, fill = Site)) + geom_point(size = 3, shape = 21) +  scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()

# Calculate the AF average of the 10 iterations of the best param
sim_AFs_best_mean <- sim_Data.melt %>%
  filter(thresh == thre, k_1 == kb1, k_2 == kb2) %>% 
  group_by(m, thresh, k_1, k_2, sim_eq) %>% reframe(AF_true_mean = mean(AF_true)) %>% 
  left_join(dplyr::select(ecovars, Site, sim_eq, ph_mean))

# Graph sim
sim_AFs_best_mean$Site <- factor(sim_AFs_best_mean$Site, levels=c("FC", "SLR", "SH", "ARA", "CBL", "PSG", "STC", "KH", "VD", "FR", "BMR", "PGP", "PL", "SBR", "PSN", "PB", "HZD", "OCT", "STR"))
pdf("output/figures/SLiM/sim_AFs_best.pdf", width = 5, height = 5)
ggplot(sim_AFs_best_mean, aes(x = AF_true_mean, y = ph_mean, fill = Site)) + geom_point(size = 3, shape = 21) + scale_fill_manual(values = mycolors) + theme_linedraw()
dev.off()