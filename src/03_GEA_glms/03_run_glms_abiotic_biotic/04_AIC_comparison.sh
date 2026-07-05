#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=AIC_comparison

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=600G 

# Submit job array
#SBATCH --array=32 #1-920%100

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 03_extract_glms.R script. 

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Determine partition to process 

# Set chunk (w) equal to slurm array
w=${SLURM_ARRAY_TASK_ID}

# Echo chunk
echo ${w}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "output_chunks" ]
then echo "Working output_chunks folder exist"; echo "Let's move on."; date
else echo "Working output_chunks folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_chunk_analysis_abiotic_biotic/output_chunks; date
fi

#--------------------------------------------------------------------------------

# Run R script
Rscript --vanilla $WORKING_FOLDER/src/03_GEA_glms/03_run_glms_abiotic_biotic/04_AIC_comparison.R "${w}"

#--------------------------------------------------------------------------------

# Say done
echo "done"