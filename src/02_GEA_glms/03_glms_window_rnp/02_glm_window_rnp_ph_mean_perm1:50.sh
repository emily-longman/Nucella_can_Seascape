#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=glms_window_rnp_perm_1:50

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=04:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=500G 

# Submit job array
#SBATCH --array=1-5 #1-995%75

# Request CPU
#SBATCH --cpus-per-task=10

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 02_glm_window_rnp_ph_mean_real.R script. 
# Note: prior to running this script, the merged glms were moved to the Working folder rather than the tmp directory.

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/wins_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. 
# Chr            nSNPs  Start   End    win_i   group
# Backbone_1553  2497  1657  101657   1        1      
# Backbone_1553  2497  51657 151657   2        1 
# ...    

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo "Window group:" ${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms/glms_window_summary

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_window_chunk_analysis_ph_mean_perm" ]
then echo "Working glms_window_chunk_analysis_ph_mean_perm folder exist"; echo "Let's move on."; date
else echo "Working glms_window_chunk_analysis_ph_mean_perm folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_window_summary/glms_window_chunk_analysis_ph_mean_perm; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript --vanilla $WORKING_FOLDER/src/02_GEA_glms/03_glms_window_rnp/02_glm_window_rnp_ph_mean_perm1:50.R "${SLURM_ARRAY_TASK_ID}"
# The --vanilla option prevents restoring or saving workspaces

#--------------------------------------------------------------------------------

# Say done
echo "done"