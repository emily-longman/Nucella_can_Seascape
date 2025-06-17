setwd("~/Documents/Seascape/biooracledata")

install.packages("gradientForest", repos="http://R-Forge.R-project.org")
####needs R version 4.4 I think? works on VACC but not on personal laptop

library(gradientForest)
library(tidyverse)
library(gameofthrones)
library(stringr)
##need to module load gdal (interpreter for raster data)
library(raster)
library(viridis) 
library(SeqArray)
library(data.table)
library(dplyr)
library(tidyr)


###load metadata (env variables extracted from biooracle)
load("indmeta.RData")

###load genotype file (creating it below if not already done)
###alternatively use AF file if population rep high enough
load("genotypefile_urchinbaypass.RData")



# Create Genotype file ----------------------------------------------------


genofile.path <- "/gpfs2/scratch/dsadler1/SeaUrchinEnvPop/GDS/Urchin_all_filtered.gds"
genofile <- seqOpen(genofile.path)

# Read significant peaks (SNPs of interest)
# Note: Ignore this step if needing all genotypes
sig_peaks <- read.table("/gpfs1/home/d/s/dsadler1/scratch/SeaUrchinEnvPop/GDS/GOwinda/subset_snps_by_peaks.txt ", col.names = c("chr", "pos"))

# Extract all SNP metadata
snps.dt <- data.table(
  variant.id = seqGetData(genofile, "variant.id"),
  chr = seqGetData(genofile, "chromosome"),
  pos = seqGetData(genofile, "position")
)

# Filter to only significant SNPs
sig_snps <- snps.dt %>%
  semi_join(sig_peaks, by = c("chr", "pos"))

# Set filter in GDS to only those SNPs
seqSetFilter(genofile, variant.id = sig_snps$variant.id)

seqGetData(genofile, "genotype") -> genotypes_i

seqGetData(genofile, "sample.id") -> sampsids_raw
gsub("/users/c/p/cpetak/sWGS/BWA_out/", "", sampsids_raw) -> sampsids_f1
gsub("_200925_A00421_0244_AHKML5DSXY_", "_", sampsids_f1) -> sampsids_f2
gsub("_L002_R1_001.rmdup.bam", "", sampsids_f2) -> sampsids_f3
gsub("18170X", "", sampsids_f3) -> sampsids_clean

colnames(genotypes_i) = sampsids_clean

genotypes_i %>% colSums() %>%
  data.frame(genotype = .) %>%
  mutate(SampleId = row.names(.)) %>%
  separate(SampleId, into = c("id","numid","Sid"), sep = "_") ->
  genotable_v1

save(genotable_v1, "genotypesallurchinGF.RData")
head(genotable_v1)

### remove NAs 
genotable_imputed <- genotable_v1
genotable_v1[is.na(genotable_v1)] <- 9

###check pop/ind matches
dim(genotable_imputed)

####probably not needed if you didnt do a stupid like me and left the ID column in...
genotable_imputed_noid<-genotable_imputed[1:4193]
head(genotable_imputed_noid)



# Running gradientforest --------------------------------------------------



###set max level i.e. maximum number of splits along any one predictor variable
###0.368  based on each tree is typically built using about 63.2% of the data, leaving 36.8% (or 0.368) out-of-bag
###prevents overfitting! (essentially removes noise from our data)
maxLevel <- log2(0.368*nrow(indmeta)/2)
maxLevel

###run gradient forest analysis 
###filters by cor threshold of 0.5 (although we have already filtered by corr matrix this is extra filtering)
####ntree can be set at 500 for testing
gfModbiooracind_new <- gradientForest(cbind(indmeta, genotable_imputed_noid), predictor.vars=colnames(indmeta),
                                  response.vars=colnames(genotable_imputed_noid), ntree=5000, 
                                  maxLevel=maxLevel, trace=T, corr.threshold=0.5)
gfModbiooracind_new



# Plotting GF outputs -----------------------------------------------------


###various plots for model metrics i.e., cumulative curves and accuracy plots
plot(gfModbiooracind_new, type="O")
plot(gfModbiooracind_new, plot.type = 'S')
plot(gfModbiooracind_new, plot.type = 'C')
plot(gfModbiooracind_new, plot.type = 'C', show.species = F)
plot(gfModbiooracind_new, plot.type = 'P')


most_important <- names(importance(gfModbiooracind_new))[1:25]

plot(gfModbiooracind_new, plot.type = "C", imp.vars = most_important,
     show.species = F, common.scale = T, cex.axis = 0.6,
     cex.lab = 0.7, line.ylab = 0.9, 
     par.args = list(mgp = c(1.5,0.5, 0),
                     mar = c(2.5, 1, 0.1, 0.5), omi = c(0,0.3, 0, 0)))

accu_imp <- as.data.frame(gfModbiooracind_new$overall.imp)
accu_imp

names(accu_imp)[1] <- "importance"
accu_imp <- tibble::rownames_to_column(accu_imp, "Variable")

accu_imp

GF_plot <- ggplot(data=accu_imp, aes(x=reorder(Variable, importance),y=(importance))) +geom_bar(stat="identity", colour="black", width=0.8)+
  theme_bw()+#scale_y_continuous(lim=c(0,15),expand = c(0, 0))+
  #Change number of colours to number of "Types" of environmental vairables -- if want to colour by type
  #scale_fill_brewer(palette="Set1")+
  scale_fill_manual(values=c("#005BFFFF","#FF0037FF"))+
  geom_hline(yintercept = 0)+
  coord_flip()+
  theme(panel.grid.major = element_blank(), legend.position=c(0.8,0.3),legend.key.size = unit(1, 'lines'),
        legend.title = element_blank(), panel.border = element_blank(),axis.ticks.y = element_blank(),
        legend.text = element_text(size = 12), axis.line.x = element_line(colour="black"), 
        panel.grid.minor = element_blank(),axis.text.x = element_text(size = 12),axis.text.y=element_text(size = 12),axis.title=element_text(size=14,color="black"),
        panel.background = element_blank())+ylab("Accuracy Importance")+ xlab("")
GF_plot


###transform present data
predOUT_present<- predict(gfModbiooracind_new, indmeta)
predOUT_present

###read in future data
futmet<-read.csv("futureindmeta.csv")
rownames(futmet)<-rownames(indmeta)

predOUT_2090 <- predict(gfModbiooracind_new, futmet)
predOUT_2090

###test it is working
(predOUT_2090[,1]-predOUT_present[,1])^2

###transform data
df_diff_squared <- data.frame(
  lapply(1:ncol(predOUT_2090), function(i) {
    (predOUT_2090[,i] - predOUT_present[,i])^2
  })
)

df_diff_squared


colnames(df_diff_squared) <- paste0("diff_sq_", colnames(predOUT_2090))
df_diff_squared

check_rsq <- gfModbiooracind_new$res
mean(check_rsq$rsq)

rownames(df_diff_squared)<-rownames(predOUT_present)

df_diff_squared <- as.data.frame(df_diff_squared) %>%
  rownames_to_column(var = "Location")
df_diff_squared



# Pool genotypes back into pops (optional) --------------------------------


####pool ind to pops (only needed for genotype level data) 
df_long_filtered <- df_diff_squared %>%
  mutate(
    Population = str_extract(Location, "^[A-Z]+")  # Extract uppercase prefix (e.g., BOD, CAP)
  ) %>%
  mutate(
    Population = factor(Population, levels = c("FOG", "CAP", "KIB", "BOD", "TER", "LOM", "SAN"))
  )

df_long_filtered <-df_long_filtered[2:9]

####


####Otherwise carry on with AF
df_long <- df_diff_squared %>%
  pivot_longer(
    cols = -Location,  # Convert all columns except Location
    names_to = "Variable",
    values_to = "Value"
  )

df_long


p1<-ggplot(df_long, aes(x = Location, y = Value, fill=Location)) +
  geom_col()+  # Use geom_col() for bar plot, or geom_point() for scatter
  facet_wrap(~ Variable, scales = "free_y") +  # Facet by bio_*
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Location", y = "Genomic offset")+
  scale_fill_got_d(option = "Tully", direction = 1)

ggsave(p1, file="p1.png")


#### change to env variables to be kept/renamed
df_filtered <- df_long %>% 
  filter(grepl("po4|ph|o2|phyc|si|sws|swd", Variable))

df_filtered <- df_long %>% 
  filter(grepl("po4|ph|o2|phyc|si|sws|swd", Variable)) %>%
  mutate(
    Variable = case_when(
      grepl("o2", Variable)   ~ "O2 (range)",
      grepl("ph", Variable)  ~ "pH (mean)",
      grepl("phyc", Variable)   ~ "Phytoplankton count (range)",
      grepl("PO4", Variable)~ "PO4 (range)",
      grepl("si", Variable)~ "Silicon content (range)",
      grepl("swd", Variable)~ "Wind direction (range)",
      grepl("sws_max", Variable)~ "Wind Speed (max)",
      grepl("sws_min", Variable)~ "Wind Speed (min)",
      grepl("sws_range", Variable)~ "Wind Speed (range)",
      
      TRUE ~ Variable  # Keep unchanged if not in the list
    )
  )
df_filtered

p2<-ggplot(df_filtered, aes(x = Location, y = Value, fill=Location)) +
  geom_col() +  # Use geom_col() for bar plot, or geom_point() for scatter
  facet_wrap(~ Variable, scales = "free_y") +  # Facet by bio_*
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Location", y = "Genomic offset")+
  scale_fill_got_d(option = "Tully", direction = 1)

ggsave(p2, file="p2.png")


df_long_filtered

####sum all columns to create overall offset values
df_diff_squared$sum_all_columns <- rowSums(df_diff_squared[2:14])


df_long_filtered

p3<-ggplot(df_diff_squared, aes(x = Location, y = sum_all_columns, fill=Location)) +
  geom_col() +  # Use geom_col() for bar plot, or geom_point() for scatter
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Location", y = "Genomic offset")+
  scale_fill_got_d(option = "Tully", direction = 1)+
  facet_wrap(~Country)

ggsave(p3, file="p3.png")


# making raster maps ------------------------------------------------------
##read in tif files generated from biooracle
envtif <- list.files(pattern = "biooracle_.*_raster.tif")
envtif <- envtif[!grepl("future", envtif, ignore.case = TRUE)]

envtiffut <- list.files(pattern = "biooracle_future.*_raster.tif")

envtif
envtiffut

###stack tif files for raster
env_stack <- stack(envtif)
env_stack_fut <- stack(envtiffut)

### Ensure proper names for layers in the stack (based on your variable names)
names(env_stack) <- sub("biooracle_(.*)_raster.tif", "\\1", envtif)
names(env_stack_fut) <- sub("biooracle_(.*)_raster.tif", "\\1", envtiffut)

### Check the layer names in the stack
print(names(env_stack))
print(names(env_stack_fut))


#####Correlation Matrix generated ######


###subset environmental variables based on correlation test (0.8 threshold)
env_stacksub <- env_stack[[
  c("o2_range", "ph_mean", "phyc_range", "po4_range", "si_max", "si_range", "swd_range", "sws_max", "sws_min", "sws_range")
]]
env_stacksub

env_stacksubfut <- env_stack_fut[[
  c("future_o2_range", "future_ph_mean", "future_phyc_range", "future_po4_range", "future_si_max", "future_si_range", "future_swd_range", "future_sws_max", "future_sws_min", "future_sws_range")
]]
env_stacksubfut

futmet2 <- futmet[
  c("o2_range", "ph_mean", "phyc_range", "po4_range", "si_max", "si_range", "swd_range", "sws_max", "sws_min", "sws_range")
]




pcaToRaster <- function(snpPreds, rast, mapCells){
  require(raster)
  
  pca <- prcomp(snpPreds, center=TRUE, scale.=FALSE)
  
  ##assigns to colors, edit as needed to maximize color contrast, etc.
  a1 <- pca$x[,1]; a2 <- pca$x[,2]; a3 <- pca$x[,3]
  r <- a1+a2; g <- -a2; b <- a3+a2-a1
  
  ##scales colors
  scalR <- (r-min(r))/(max(r)-min(r))*255
  scalG <- (g-min(g))/(max(g)-min(g))*255
  scalB <- (b-min(b))/(max(b)-min(b))*255
  
  ##assigns color to raster
  rast1 <- rast2 <- rast3 <- rast
  rast1[mapCells] <- scalR
  rast2[mapCells] <- scalG
  rast3[mapCells] <- scalB
  ##stacks color rasters
  outRast <- stack(rast1, rast2, rast3)
  return(outRast)
}

# Function to map difference between spatial genetic predictions
# predMap1 = dataframe of transformed variables from gf or gdm model for first set of SNPs
# predMap2 = dataframe of transformed variables from gf or gdm model for second set of SNPs
# rast = a raster mask to which Procrustes residuals are to be mapped
# mapCells = cell IDs to which Procrustes residuals values should be assigned
RGBdiffMap <- function(predMap1, predMap2, rast, mapCells){
  require(vegan)
  PCA1 <- prcomp(predMap1, center=TRUE, scale.=FALSE)
  PCA2 <- prcomp(predMap2, center=TRUE, scale.=FALSE)
  diffProcrust <- procrustes(PCA1, PCA2, scale=TRUE, symmetrical=FALSE)
  residMap <- residuals(diffProcrust)
  rast[mapCells] <- residMap
  return(list(max(residMap), rast))
}

# OK, on to mapping. Script assumes:
# (1) a dataframe named env_trns containing extracted raster data (w/ cell IDs)
# and env. variables used in the models & with columns as follows: cell, bio1, bio2, etc.
#
# (2) a raster mask of the study region to which the RGB data will be written

# transform env using gf models, see ?predict.gradientForest

env_df <- as.data.frame(env_stack, xy = TRUE, cell = TRUE)
env_df_fut <- as.data.frame(env_stacksubfut, xy = TRUE, cell = TRUE)

env_df$cell <- cellFromXY(env_stack, env_df[, c("x", "y")])
env_df_fut$cell <- cellFromXY(env_stacksubfut, env_df_fut[, c("x", "y")])


# Reorder columns to place "cell" first
env_df <- env_df %>% dplyr::select(cell, everything())
env_df
env_dfsub<-env_df %>% 
  dplyr::select("o2_range", "ph_mean", "phyc_range", "po4_range", "si_max", "si_range", "swd_range", "sws_max", "sws_min", "sws_range")


env_df <- na.omit(env_df)
env_df

env_df_fut <- na.ostrmit(env_df_fut)


colnames(env_df_fut)<-colnames(env_df)
env_df_fut
env_dfsub <- na.omit(env_dfsub)
env_dfsub


env_dfsubfut<-env_df_fut %>% 
  dplyr::select("future_o2_range", "future_ph_mean", "future_phyc_range", "future_po4_range", "future_si_max", "future_si_range", "future_swd_range", "future_sws_max", "future_sws_min", "future_sws_range")

colnames(env_dfsubfut)<-colnames(env_dfsub)
env_dfsubfut


##transform based on full range
predRefind <- predict(gfModbiooracind_new, env_dfsub) # remove cell column before transforming
predRefindfut <- predict(gfModbiooracind_new, env_dfsubfut) # remove cell column before transforming


###load blank raster file
load("study_raster.RData")
study_raster
mask <- study_raster
mask

###or create as such: 
latitude_range <- c(30, 45)
longitude_range <- c(-125, -112)

study_extent <- extent(longitude_range[1], longitude_range[2], latitude_range[1], latitude_range[2])

raster_resolution <- 0.05

study_raster <- raster(study_extent, res = raster_resolution)

crs(study_raster) <- CRS("+proj=longlat +datum=WGS84")
values(study_raster) <- NA
plot(study_raster, main = "Raster Template for California and Oregon")

###present data map
refRGBmap <- pcaToRaster(predRefind, mask, env_df$cell)
plotRGB(refRGBmap)
writeRaster(refRGBmap, "/.../refSNPs_map.tif", format="GTiff", overwrite=TRUE)

#futuredata
GI5RGBmap <- pcaToRaster(predRefindfut, mask, env_df_fut$cell)
plotRGB(GI5RGBmap)
writeRaster(refRGBmap, "/.../GI5SNPs_map.tif", format="GTiff", overwrite=TRUE)

# Difference between maps (future and present) 
diffGI5 <- RGBdiffMap(predRefind, predRefindfut, rast=mask, mapCells=env_df$cell)
plot(diffGI5[[2]])
writeRaster(diffGI5[[2]], "/.../diffRef_GI5.tif", format="GTiff", overwrite=TRUE)



# calculate euclidean distance between current and future genetic spaces  
genOffsetGI5 <- sqrt((predRefindfut[,1]-predRefind[,1])^2+(predRefindfut[,2]-predRefind[,2])^2
                     +(predRefindfut[,3]-predRefind[,3])^2+(predRefindfut[,4]-predRefind[,4])^2
                     +(predRefindfut[,5]-predRefind[,5])^2+(predRefindfut[,6]-predRefind[,6])^2
                     +(predRefindfut[,7]-predRefind[,7])^2)


# assign values to raster - can be tricky if current/future climate
# rasters are not identical in terms of # cells, extent, etc.
mask
mask[env_df_fut$cell] <- genOffsetGI5
plot(mask)

plot(diffGI5[[2]])
pdf("map.pdf", width = 8, height = 6)
plot(diffGI5[[2]])
dev.off()

pdf("mask.pdf", width = 8, height = 6)
plot(mask)
dev.off()


# Convert raster to a dataframe for ggplot
mask_df <- as.data.frame(rasterToPoints(mask))
colnames(mask_df) <- c("Longitude", "Latitude", "GenOffset")

# Create a better plot with ggplot2
p1<-ggplot(mask_df, aes(x = Longitude, y = Latitude, fill = GenOffset)) +
  geom_raster() +
  scale_fill_viridis(option = "magma", na.value = "transparent", name = "Genomic Offset") + 
  theme_minimal() +
  coord_fixed() +  # Keeps aspect ratio
  labs(title = "Genomic Offset Map", 
       x = "Longitude", 
       y = "Latitude") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )



mask_highres <- disaggregate(mask, fact = 2)  # Increase resolution 2x
mask_df <- as.data.frame(rasterToPoints(mask_highres))
mask_df
p1<-ggplot(mask_df, aes(x = x, y = y, fill = layer)) +
  geom_tile() +  # Instead of geom_raster()
  scale_fill_viridis(option = "magma", na.value = "transparent", name = "Genomic Offset") + 
  theme_minimal() +
  coord_fixed() +
  labs(title = "Genomic Offset Map", x = "Longitude", y = "Latitude")

ggsave(p1, file="mapkelpGF.png")




# extra stuff ignore for now... -------------------------------------------


setwd("~/Documents/GenomicOffset/")
climGF2<-read.csv(file = "MetadataNASACOP.csv")
climGF2<-climGF2[,-1]
rownames(climGF2)<-climGF2$id
climGF2<-climGF2[,-1]
climGF2 <- climGF2[, !names(climGF2) %in% "PRECTOTCORR_min"]
climGF2

bioclimGFfuture<-read.csv("bioclimfuture.csv", header=T)
bioclimGFfuture
rownames(bioclimGFfuture)<-c("FOG", "CAP", "KIB", "BOD", "TER", "LOM", "SAN")
bioclimGFfuture<-bioclimGFfuture[,-1]
bioclimGFfuture<-bioclimGFfuture[,-1]
bioclimGFfuture
write.csv(bioclimGFfuture, file="bioclimfutureformatted.csv")




bioclimGF<-read.csv("bioclimpresent.csv", header=T)
bioclimGF
rownames(bioclimGF)<-c("FOG", "CAP", "KIB", "BOD", "TER", "LOM", "SAN")
bioclimGF<-bioclimGF[,-1]
bioclimGF<-bioclimGF[,-1]
bioclimGF
write.csv(bioclimGF, file="bioclimpresentformatted.csv")

biooracGF<-read.csv("allpresentdatasub.csv", header=T)
dim(biooracGF)
biooracGF
biooracGF<-biooracGF %>% filter(time=="2010-01-01T00:00:00Z")
biooracGF
dim(biooracGF)
bio<-biooracGF[, -c(1,16)]
bio
cor_site_values<-cor(bio, method = "pearson")
cor_site_values
write.csv(cor_site_values,file="Correlation_site_values.csv")




# Identify highly correlated pairs (e.g., |correlation| > 0.8)
high_corr <- findCorrelation(cor_site_values, cutoff = 0.8)

high_corr
# Remove highly correlated variables
cor_site_values_filter <- cor_site_values[, -high_corr]
cor_site_values_filter


env1_cor <- cor(bio,method = "pearson")
env1_cor
library(corrplot)

corrplot(cor_site_values)
str(bio)

model <- lm(thetao_max + thetao_min + thetao_range + thetao_mean + chl_mean + no3_mean +
              o2_mean + ph_max + ph_min + ph_range + ph_mean + phyc_mean + po4_mean + so_mean ~ ., 
            data = bio)
model

vif(model_max)

model_max <- lm(thetao_max ~ thetao_min + thetao_range + thetao_mean + chl_mean + 
                  no3_mean + o2_mean + ph_max + ph_min + ph_range + ph_mean + 
                  phyc_mean + po4_mean + so_mean, data = bio)
model_max


cor_matrix <- cor(bio[, c("thetao_min", "thetao_range", "thetao_mean", "chl_mean", 
                          "no3_mean", "o2_mean", "ph_max", "ph_min", "ph_range", 
                          "ph_mean", "phyc_mean", "po4_mean", "so_mean")])