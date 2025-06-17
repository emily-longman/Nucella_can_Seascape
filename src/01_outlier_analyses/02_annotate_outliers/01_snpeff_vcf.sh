#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=snpeff

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run snpeff on the vcf file. 

# Load modules
module load snpeff/5.2c

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# VCF is the input vcf file
VCF=$WORKING_FOLDER/data/processed/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter_minQ60_maxmissing1.0.recode.vcf

# This is the location where the snpeff directory was built.
DATA_DIR=/netfiles/pespenilab_share/Nucella/processed

# This is the snpeff config file.
PARAM=/netfiles/pespenilab_share/Nucella/processed/N.canaliculata_snpeff_March_2025/snpEff.config

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/outlier_analyses

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "snpeff" ]
then echo "Working snpeff folder exist"; echo "Let's move on."; date
else echo "Working snpeff folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/snpeff; date
fi

#--------------------------------------------------------------------------------

# Run SNPeff on the vcf file

# Change directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/snpeff

# Run SNPeff
snpeff -c $PARAM -dataDir $DATA_DIR N.canaliculata_snpeff_March_2025 $VCF > $WORKING_FOLDER/data/processed/outlier_analyses/snpeff/N.canaliculata_pops_SNPs_annotate.vcf