#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=model_enrichment

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=02:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=300G 

# Submit job array
#SBATCH --array=1-9 #10-27

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 02_model_enrichment.R. 

# Load modules 
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/Seascape_vars_names.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. 
# Env variable name   
# thetao_max
#...    

#--------------------------------------------------------------------------------

# Determine the environmental variable to process

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $guide_file`

# State environmental variable
echo ${i}

#--------------------------------------------------------------------------------

# Run R script

Rscript --vanilla $WORKING_FOLDER/src/02_GEA_glms/02_model_enrichment/03_model_enrichment.R "${i}"

#--------------------------------------------------------------------------------

# Say done
echo "done"