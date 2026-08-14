# Generate PCA of all SNPs

# Clear memory
rm(list=ls()) 

# ================================================================================== #

# Note: to run script move gwas output to results dir 

# Set path as main Github repo
install.packages(c('rprojroot'))
library(rprojroot)

# List all files and directories below the root
dir(find_root(has_file("README.md")))

# Set relative path of results directory from root
dir(find_root_file("data", "processed",  criterion = has_file("README.md")))
data_path_from_root <- find_root_file("data", "processed", criterion = has_file("README.md"))
# List files in this folder to make sure you're in the right spot.
list.files(data_path_from_root)

# Set working directory as path from root
setwd(data_path_from_root)

# ================================================================================== #

# Load packages
install.packages(c('adegenet', 'graph4lg'))
require(adegenet)
require(graph4lg)

# ================================================================================== #

# Load Data: 

# Load SNP data (rows are populations and each column is a lcous)
data <- read.table("pop_structure/genotype_table/N.can.pop.genotypes.txt", header = F, sep = "\t", stringsAsFactors=F)
# Get population names
pops <- read.table("pop_structure/genotype_table/Nucella_pops.list", header=F)

# Reformat data

# Rename rownames as SNP names, then remove first column which contained the rownames
row.names(data) = data$V1
data=data[,-1]

dim(data) # 265  19
str(data)

# ================================================================================== #

# How many alleles per genotype 

datcount=t(apply(data,1,function(x){xx=as.numeric(strsplit(paste(x,collapse="/"),"/")[[1]]);sum(!is.na(xx) & !duplicated(xx),na.rm=T);}))
datcount=t(datcount)
head(datcount)
table(datcount)



# ================================================================================== #


# Graph pca for 1 M random SNPs

pca=prcomp(t(data), scale=T)
