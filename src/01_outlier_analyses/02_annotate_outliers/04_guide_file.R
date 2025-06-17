#### Create guide file for array via an interactive session.

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# PROTEIN_FILE is the N. canaliculata protein file generated from the gff
PROTEIN_FILE=/netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/protein.fa

#--------------------------------------------------------------------------------

# Currently the protein file is interleaven (i.e., the DNA sequence is in chunks of 80bps a line).
# For easy of running a loop and array, make it not interleaven

# Load modules 
seqtk=/gpfs1/home/e/l/elongman/software/seqtk/seqtk

# Use seqtk to make protein file not interleaven
$seqtk seq $PROTEIN_FILE > $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/Nucella_canaliculata_protein.fa

#--------------------------------------------------------------------------------

# Get list of protein names
 
# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology

# Make guide files directory
mkdir guide_files

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/guide_files

# Get all of the scaffold/contig names 
grep ">" $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/Nucella_canaliculata_protein.fa > protein_file_names.txt

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------

# Modify in R
module load R/4.4.1  
R

# Load packages
install.packages(c('data.table', 'tidyverse', 'groupdata2'))
library(data.table)
library(tidyverse)
library(groupdata2)

#--------------------------------------------------------------------------------

# Read in the data 
data <- fread("protein_file_names.txt", header = F)
dim(data) # 40634  4

# Break up the data into 3 runs - will run a bash script for each run
data_run1 <- data[1:14000,]
data_run2 <- data[14001:28000,]
data_run3 <- data[28001:40634,]

# Break up the protein names up into chunks
group(data_run1, n=15, method = "greedy") -> protein_file_names_run_1
group(data_run2, n=15, method = "greedy") -> protein_file_names_run_2
group(data_run3, n=15, method = "greedy") -> protein_file_names_run_3

# Write guide file
write.table(protein_file_names_run_1, "protein_file_names_array_run_1.txt", col.names = F, row.names = F, quote = F)
write.table(protein_file_names_run_2, "protein_file_names_array_run_2.txt", col.names = F, row.names = F, quote = F)
write.table(protein_file_names_run_3, "protein_file_names_array_run_3.txt", col.names = F, row.names = F, quote = F)

# Note each guide file has dimensions:  14000, 5 (934 partitions with 15 protein names in each)

q()