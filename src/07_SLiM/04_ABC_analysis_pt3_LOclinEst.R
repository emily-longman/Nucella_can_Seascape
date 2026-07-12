## Nucella ABC analysis
## Part 3 - ABC estimation

library(tidyverse)
library(data.table)
library(magrittr)
library(reshape2)
library(gmodels)
library(poolfstat)
library(foreach)
library(lme4)
library(abc)

###
#load
real_All <- get(load("/gpfs2/scratch/jcnunez/Nucella_Sims_EcoLoad/real.data.Rdata"))
sim_All <- get(load("/gpfs2/scratch/jcnunez/Nucella_Sims_EcoLoad/sim.results.Rdata"))
### estimate means

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
### .. EXCLUDE m and N as those are invariant

###
abc_fit <- abc(
  target = real_data,
  param = sim_parameters,
  sumstat = simulated_data,
  tol = 0.1,
  method = "loclinear"
)
best <- which.max(abc_fit$weights)
abc_fit$unadj.values[best, ]

post <- as.data.frame(abc_fit$adj.values)
post_long <- pivot_longer(
  post,
  cols = everything(),
  names_to = "parameter",
  values_to = "value"
)
post_long

ggplot(post_long, aes(x = value)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  facet_wrap(parameter~. , scales = "free", ncol = 1) +
  theme_bw() -> posteriors
ggsave(posteriors, file = "posteriors.pdf")

### Now report the simulation that better match the data
post_long %>%
  group_by(parameter) %>%
  summarise(mean.par = mean(value, na.rm = T))


best <- which.max(abc_fit$weights)
abc_fit$unadj.values[best, ]

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

kb1 = abc_fit$unadj.values[best,"k_1"]
kb2 = abc_fit$unadj.values[best,"k_2"]
thre = abc_fit$unadj.values[best,"thresh"]

f <- function(x, z, k) {
  2/(1+exp(-(x - z)/k))-1
}

data.frame(env=env, s=f(env, thre, kb1)) %>%
  ggplot(aes(
    x=env, y=s
  )) + geom_line() ->
  sel_curve

ggsave(sel_curve, file = "sel_curve.pdf")

### Compare to real data
### Compare to real data
### Compare to real data
### Compare to real data

###load ecovars
ecovars <- fread("/gpfs2/scratch/jcnunez/Nucella_Sims_EcoLoad/Nucella_ph_shellt.txt")
names(ecovars)[2] = "Site"
ecovars %<>% mutate(sim_eq = paste("p", 0:18, sep =""))

##load top SNP
phafs <- fread("/gpfs2/scratch/jcnunez/Nucella_Sims_EcoLoad/afs.ph.g27343.BF.POD.csv") %>%
  mutate(nsnails = 20)
phafs %>% 
  filter(SNP_id == "ntLink_3821_1595") %>%
  left_join(dplyr::select(ecovars, Site, ph_mean))->
  topsnp

####
sim_Data <- get(load("/netfiles/nunezlab/Shared_Resources/in_transit/longman/ph_results.Rdata"))
names(sim_Data)[1:19]= paste("p", 0:18, sep ="")
sim_Data %>%
  reshape2::melt(id = c("repId","m","thresh","k_1","k_2",
                        "N","state","sim.cycle"),
                 variable.name = "sim_eq",
                 value.name = "AF_true") -> 
  sim_Data.melt


###recall...
kb1 = abc_fit$unadj.values[best,"k_1"]
kb2 = abc_fit$unadj.values[best,"k_2"]
thre = abc_fit$unadj.values[best,"thresh"]

w <- abc_fit$weights / sum(abc_fit$weights)
colSums(
  abc_fit$unadj.values * w
)



sim_Data.melt %>% 
  filter( 
          thresh == thre,
          k_1 == kb1,
          k_2 ==kb2 
          ) %>%
  left_join(dplyr::select(ecovars, sim_eq, ph_mean)) %>%
  ggplot(aes(
    y=ph_mean, x=1-AF_true
  )) + geom_point(size = 2.0) ->
  sim_AFs

ggsave(sim_AFs, file = "sim_AFs.pdf")

