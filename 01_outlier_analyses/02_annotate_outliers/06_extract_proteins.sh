#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=extract_blast_proteins

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=15G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will extract the protein data for the SNPs of interest.

# Load modules
seqtk=/gpfs1/home/e/l/elongman/software/seqtk/seqtk
ncbi=/gpfs1/home/e/l/elongman/software/ncbi-blast-2.16.0+/bin

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/outlier_analyses

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "annotation" ]
then echo "Working annotation folder exist"; echo "Let's move on."; date
else echo "Working annotation folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/annotation; date
fi

#--------------------------------------------------------------------------------

# Move a copy of the protein file and make it not interleaved

# Use seqtk to make the protein file not interleaved
#$seqtk seq /netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/protein.fa \
#> $WORKING_FOLDER/data/processed/outlier_analyses/protein.fa

#--------------------------------------------------------------------------------

# SNP list (open csv then resave as txt to remove "")
SNP=$WORKING_FOLDER/data/processed/outlier_analyses/snpeff/N.canaliculata_SNPs_bonferroni_annotated_mod.txt

# Extract unique protein names
awk -F "\t" '{print $11}' $SNP | sed '1d' | uniq > $WORKING_FOLDER/data/processed/outlier_analyses/annotation/SNP.protein.names.txt

#--------------------------------------------------------------------------------

# Extract and blast proteins of outlier SNPs
cat $WORKING_FOLDER/data/processed/outlier_analyses/annotation/SNP.protein.names.txt | \
while read prot
do echo ${prot}

grep -EA 1 "${prot}" $WORKING_FOLDER/data/processed/outlier_analyses/protein.fa > $WORKING_FOLDER/data/processed/outlier_analyses/annotation/${prot}.fa

# Use the blastp command to compare the Nucella protein file with the uniprot database
$ncbi/blastp -query $WORKING_FOLDER/data/processed/outlier_analyses/annotation/${prot}.fa \
-db $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/uniprot/uniref90 \
-out $WORKING_FOLDER/data/processed/outlier_analyses/gene_ontology/annotation \
-outfmt 6 \
-evalue 1e-5 \
-max_target_seqs 5

# query <File_In>: Input file name
# db: BLAST database name
# outfmt 6: alignment view options: Tabular
# evalue: Expectation value (E) threshold for saving hits. Default = 10 
# max_target_seqs: Maximum number of aligned sequences to keep 

done
