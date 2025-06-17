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

# Submit job array
#SBATCH --array=999 #0-999%50

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 02_glms_bio-oracle.R script. This will run glms on the 

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/05_GEA/02_glms

#--------------------------------------------------------------------------------

# Input files
chunk=1
array=$((${SLURM_ARRAY_TASK_ID}+1))

echo "I am running chunk:" ${chunk} "and array:" ${array}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_window_analysis" ]
then echo "Working glms_window_analysis folder exist"; echo "Let's move on."; date
else echo "Working glms_window_analysis folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_window_analysis; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/GEA/glms/glms_window_analysis

if [ -d "GLM_100perm_Bio-Oracle_chunk_${chunk}" ]
then echo "Working GLM_100perm_Bio-Oracle_chunk_${chunk} folder exist"; echo "Let's move on."; date
else echo "Working GLM_100perm_Bio-Oracle_chunk_${chunk} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_window_analysis/GLM_100perm_Bio-Oracle_chunk_${chunk}; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript $SCRIPT_FOLDER/02_glms_bio-oracle.R "$chunk" "$array"