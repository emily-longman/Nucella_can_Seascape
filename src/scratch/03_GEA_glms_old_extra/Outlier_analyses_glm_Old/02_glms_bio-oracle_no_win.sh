#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=glms_bio-oracle

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=1:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Request CPU
#SBATCH --cpus-per-task=8

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 02_glms_bio-oracle.R script. 
# This will run glms and look at the association of the allele frequencies of the outlier SNPs and the Bio-oracle environmental data.

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/02_GEA/02_glms

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_no_window" ]
then echo "Working glms_no_window folder exist"; echo "Let's move on."; date
else echo "Working glms_no_window folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_no_window; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/GEA/glms/glms_no_window

if [ -d "GLM_100perm_Bio-Oracle_chunk_${chunk}" ]
then echo "Working GLM_100perm_Bio-Oracle_chunk_${chunk} folder exist"; echo "Let's move on."; date
else echo "Working GLM_100perm_Bio-Oracle_chunk_${chunk} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_no_window; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript $SCRIPT_FOLDER/02_glms_bio-oracle_no_win.R