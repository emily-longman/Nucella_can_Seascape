#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=merge_glms_biotic_76:100

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=8:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=900G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 03_merge_glms.R script. 

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/02_GEA/02_glms/01_run_glms_biotic

#--------------------------------------------------------------------------------

# Run R script

#Rscript --vanilla $SCRIPT_FOLDER/02_merge_glms_real.R
#Rscript --vanilla $SCRIPT_FOLDER/02_merge_glms_perm_1:25.R
#Rscript --vanilla $SCRIPT_FOLDER/02_merge_glms_perm_26:50.R
#Rscript --vanilla $SCRIPT_FOLDER/02_merge_glms_perm_51:75.R
Rscript --vanilla $SCRIPT_FOLDER/02_merge_glms_perm_76:100.R
# The --vanilla option prevents restoring or saving workspaces

#--------------------------------------------------------------------------------

# Say done
echo "done"