#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=variant_calling_bcftools_cat

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will cat together the chunked vcf files generated in step 08 from freebayes.

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Combine together the vcf files from all of the partitions

bcftools concat \
$WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.*.vcf.gz \
-Oz -o $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz