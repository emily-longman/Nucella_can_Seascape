#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=LFMM_biotic_MARINe_perm1:25

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=200G 

# Request CPU
#SBATCH --cpus-per-task=20

# Submit job array
#SBATCH --array=1-14

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 03_LFMM_biotic_MARINe_perm.R script. 

# Load modules 
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Guide file 
ECOL_VARIABLES=$WORKING_FOLDER/guide_files/Seascape_vars_names_biotic_MARINe.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. 
# Env variable name   
# thetao_max
#...    

#--------------------------------------------------------------------------------

# Determine the ecological variable to process

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $ECOL_VARIABLES`

# State ecological variable
echo ${i}

#--------------------------------------------------------------------------------

# Run R script
Rscript --vanilla $WORKING_FOLDER/src/03_LFMM/03_LFMM_biotic_MARINe_perm1:25.R "${i}"

#--------------------------------------------------------------------------------

# Say done
echo "done"