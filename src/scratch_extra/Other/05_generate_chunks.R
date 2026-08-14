# Generate 1000 chunks for baypass

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
install.packages(c('data.table', 'tidyverse'))
library(data.table)
library(tidyverse)

# ================================================================================== #

# Set seed for reproducibility
set.seed(123)

# Path to files
#genotype_file <- "data/processed/fastq_to_vcf/genotype_table/N.canaliculata.genotypes.txt" 
genobaypass_file <- "data/processed/GEA/baypass/genobaypass" 
#snpdet_file <- "data/processed/GEA/baypass/snpdet" 
output_dir <- "data/processed/GEA/baypass/baypass_chunks"

# Ensure output directory exists
if (!dir.exists(output_dir)) dir.create(output_dir)

# ================================================================================== #

# Load the genotype file
#geno_data <- fread(genotype_file, header = FALSE)
#str(geno_data)
genobaypass <- fread(genobaypass_file, sep=" ", header=F)
str(genobaypass)

# Assign row numbers as SNP IDs
snp_ids <- seq_len(nrow(genobaypass))  
# Get total number of SNPs (8,277,206)
num_snps <- nrow(genobaypass) 

# ================================================================================== #

# Load the snp info matrix from poolfstat2baypass
# Note: we need to use this SNP list since we used poolfstat to do more filtering
#snpdet <- read.table(snpdet_file, sep="", header=F)
#str(snpdet)
#dim(snpdet)

# Create new snpdet column that has both chr and position separated by a space
#snpdet$CHR_POS <- paste(snpdet$V1, snpdet$V2, sep=" ")

# ================================================================================== #

# Subset the genotype file based on the SNPs in snpdet
#geno_data_sub <- dplyr::semi_join(geno_data, snpdet, by = c('V1' = 'CHR_POS'))

# Check to make sure the number of SNPs match
#dim(geno_data_sub)[1] #8277160

# Assign row numbers as SNP IDs
#snp_ids <- seq_len(nrow(geno_data_sub))  
# Get total number of SNPs
#num_snps <- nrow(geno_data_sub) 

# ================================================================================== #

# Set number of chunks
num_chunks <- 1000

# Shuffle SNP indices
random_indices <- sample(1:num_snps, num_snps)

# Split into 1000 chunks
chunk_size <- ceiling(num_snps / num_chunks) #8,278
chunks <- split(random_indices, ceiling(seq_along(random_indices) / chunk_size))

# Save each chunk with SNP IDs separately
for (i in seq_along(chunks)) {
  chunk_indices <- chunks[[i]]
  chunk_data <- genobaypass[chunk_indices, ]  # Extract shuffled genotype data
  chunk_ids <- snp_ids[chunk_indices]  # Corresponding row numbers as SNP IDs
  
  # Double-check column count before writing
  if (ncol(chunk_data) != 38) {
    stop(paste("Error: Chunk", i, "does not have 38 columns!"))
  }
  
  # Save genotype data (BayPass input)
  fwrite(chunk_data, 
         file = sprintf("%s/genobaypass_chunk_%04d.txt", output_dir, i), 
         row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Save SNP IDs for reordering later
  fwrite(data.table(Index = chunk_indices, SNP_ID = chunk_ids), 
         file = sprintf("%s/metadata_chunk_%04d.txt", output_dir, i), 
         row.names = FALSE, col.names = TRUE, sep = "\t")
}