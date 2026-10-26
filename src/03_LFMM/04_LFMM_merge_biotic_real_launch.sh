#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=LFMM_merge_real_launch_biotic

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=01:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 04_LFMM_merge_launch.R script. 

# Load modules 
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Run R script
Rscript --vanilla $WORKING_FOLDER/src/03_LFMM/04_LFMM_merge_abiotic_real.R

#--------------------------------------------------------------------------------

# Say done
echo "done"