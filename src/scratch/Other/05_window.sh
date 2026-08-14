#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_windows

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=3:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=15G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# Generate windows for baypass analyses.

#Load modules 
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# SCRIPT_FOLDER is the folder where the scripts are.
SCRIPT_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/src

#--------------------------------------------------------------------------------

# Define input files
XtX_file=$WORKING_FOLDER/GEA/baypass/xtx/NC_baypass_core_summary_pi_xtx.out
SNP_meta_file=$WORKING_FOLDER/GEA/baypass/snpdet

# Run R script
Rscript $SCRIPT_FOLDER/03_GEA/02_baypass/05_window.R "$XtX_file" "$SNP_meta_file"
