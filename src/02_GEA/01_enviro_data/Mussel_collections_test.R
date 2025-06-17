# ------------------------------------------------------------------------------
# Mussel collections estimate:  
# Determine the number of mussels and size dist needed to get a representative sample
# Emily K. Longman
# ------------------------------------------------------------------------------

# Load libraries
library(ggplot2)
library(dplyr)
library(tidyverse)
library(here)

# Load data --------------------------------------------------------------------
mussels <- read.csv(here::here("data/raw/enviro/M.cali.thk/Mcal_shell_thickness_2019.csv"))
mussels$Thk.length <- mussels$Ave.thk.1.3/ mussels$Shell.Length
str(mussels)

# Look at all data -------------------------------------------------------------

# Graph all mussels
ggplot(mussels, aes(Shell.Length, Ave.thk.1.3)) + geom_point()

#Linear regression
mussel.mod <- lm(Ave.thk.1.3 ~ Shell.Length, mussels)
plot(mussel.mod)
summary(mussel.mod)

#calc thickness/ length ratio
mean(mussels$Thk.length) #0.02117015

#Just Bodega Mussels -----------------------------------------------------------
mussels.BMR <- mussels[which(mussels$Site.Code == "BH"),]

length(mussels.BMR) #101 mussels

#calc thickness/ length ratio
mean(mussels.BMR$Thk.length) #0.02124721

ggplot(mussels.BMR, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.BMR <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.BMR)
summary(mussel.mod.BMR)

#Sample just 50 BMR mussels across all sizes 100 times 
mussels.BMR.sample = NULL
for (i in 1:100) {
  s = sample(mussels.BMR$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.BMR.sample<- c(mussels.BMR.sample, m)
  }
hist(mussels.BMR.sample)
abline(v = 0.02124721, col= 2)

range(mussels.BMR.sample)


#Sample 50 BMR mussels that are less than 120mm 100 times 
mussels.BMR.less.120 <- mussels.BMR[which(mussels.BMR$Shell.Length < 120),]

mussels.BMR.less.120.sample = NULL
for (i in 1:100) {
  s = sample(mussels.BMR.less.120$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.BMR.less.120.sample<- c(mussels.BMR.less.120.sample, m)
  }
hist(mussels.BMR.less.120.sample)
abline(v = 0.02124721, col= 2)

range(mussels.BMR.less.120.sample)


# Just Strawberry Hill Mussels -------------------------------------------------
mussels.SH <- mussels[which(mussels$Site.Code == "SH"),]

dim(mussels.SH) #108 mussels
mean(mussels.SH$Ave.thk.1.3/ mussels.SH$Shell.Length) #0.02396629

ggplot(mussels.SH, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.SH <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.SH)
plot(mussel.mod.SH)
summary(mussel.mod.SH)

#Sample just 50 SH mussels across all sizes 100 times 
mussels.SH.sample = NULL
for (i in 1:100) {
  s = sample(mussels.SH$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.SH.sample<- c(mussels.SH.sample, m)
  }
hist(mussels.SH.sample)
abline(v = 0.02396629, col= 2)

range(mussels.SH.sample) #0.02275464 0.02591909

#Sample 50 SH mussels that are less than 120mm and greater than 40mm 100 times 
mussels.SH.less.120 <- mussels.SH[which(mussels.SH$Shell.Length < 120 & mussels.SH$Shell.Length > 40),]
mussels.SH.less.120.sample = NULL
for (i in 1:100) {
  s = sample(mussels.SH.less.120$Thk.length, 50, replace=TRUE)
  m = mean(s)
  mussels.SH.less.120.sample<- c(mussels.SH.less.120.sample, m)
  }
hist(mussels.SH.less.120.sample)
abline(v = 0.02396629, col= 2)

range(mussels.SH.less.120.sample)

# Just Fogarty Creek Mussels ---------------------------------------------------
mussels.FC <- mussels[which(mussels$Site.Code == "FC"),]

dim(mussels.FC) #108 mussels
mean(mussels.FC$Ave.thk.1.3/ mussels.FC$Shell.Length) #0.01994755

ggplot(mussels.FC, 
       aes(Shell.Length, Ave.thk.1.3)) + geom_point()

mussel.mod.FC <- lm(Ave.thk.1.3 ~ Shell.Length, mussels.FC)
plot(mussel.mod.FC)
summary(mussel.mod.FC)

#Sample just 50 SH mussels across all sizes 100 times 
mussels.FC.sample = NULL
for (i in 1:100) {
  s = sample(mussels.FC$Thk.length, 55, replace=TRUE)
  m = mean(s)
  mussels.FC.sample<- c(mussels.FC.sample, m)
  }
hist(mussels.FC.sample)
abline(v = 0.01994755, col= 2)

range(mussels.FC.sample) 

#Sample just 50 SH mussels that are <120 and greater than 40mm all sizes 100 times 
mussels.FC.less.120.great.40 <- mussels.FC[which(mussels.FC$Shell.Length < 120 
                                                 & mussels.FC$Shell.Length > 40),]
mussels.FC.less.120.great.40.sample = NULL
for (i in 1:100) {
  s = sample(mussels.FC.less.120.great.40$Thk.length, 55, replace=TRUE)
  m = mean(s)
  mussels.FC.less.120.great.40.sample<- c(mussels.FC.less.120.great.40.sample, m)
  }
hist(mussels.FC.less.120.great.40.sample)
abline(v = 0.01994755, col= 2)

range(mussels.FC.less.120.great.40.sample) 



