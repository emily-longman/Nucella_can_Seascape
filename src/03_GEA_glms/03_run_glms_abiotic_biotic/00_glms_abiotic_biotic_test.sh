#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=glms_abiotic_biotic_test

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=15:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Request CPU
#SBATCH --cpus-per-task=10


# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 01_glms_Mcali.R script.  
# This will run glms and look at the association of the allele frequencies of the outlier SNPs and the Marine biotic data.

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/03_GEA_glms/03_run_glms_abiotic_biotic

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_chunk_analysis_abiotic_biotic" ]
then echo "Working glms_chunk_analysis_abiotic_biotic folder exist"; echo "Let's move on."; date
else echo "Working glms_chunk_analysis_abiotic_biotic folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript $SCRIPT_FOLDER/00_glms_abiotic_biotic_perm1:5.R

#--------------------------------------------------------------------------------

# Say done
echo "done"