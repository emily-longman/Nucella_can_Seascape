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
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Request CPU
#SBATCH --cpus-per-task=5

# Submit job array
#SBATCH --array=1-500%50

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

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
SCRIPT_FOLDER=$WORKING_FOLDER/src/02_GEA_glms/01_run_glms_abiotic

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/scaffold_names_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 18919, 2; 996 partitions each with 19 scaffold names)
# Scaffold name         # Partition/array
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# List scaffold names
cat scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_chunk_analysis_abiotic" ]
then echo "Working glms_chunk_analysis_abiotic folder exist"; echo "Let's move on."; date
else echo "Working glms_chunk_analysis_abiotic folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_chunk_analysis_abiotic; date
fi

#--------------------------------------------------------------------------------

# Run R script

Rscript $SCRIPT_FOLDER/02_glms_bio-oracle.R "${SLURM_ARRAY_TASK_ID}"

#--------------------------------------------------------------------------------

rm scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# Say done
echo "done"