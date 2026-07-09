#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ph_slim

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Submit job array
#SBATCH --array=1-2

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying slim script. 

# Load modules 
module load slim/5.0

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

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
cd $WORKING_FOLDER/data/processed

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "SLiM" ]
then echo "Working SLiM folder exist"; echo "Let's move on."; date
else echo "Working SLiM folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "results" ]
then echo "Working results folder exist"; echo "Let's move on."; date
else echo "Working results folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/results; date
fi

#--------------------------------------------------------------------------------

repid=${SLURM_ARRAY_TASK_ID}
root=$WORKING_FOLDER/data/processed/SLiM/results

# Guide file 
#guide_file=$WORKING_FOLDER/guide_files/scaffold_names_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
# Scaffold name         # Partition/array
# ....

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/results

# Run slim script
slim 	
-d "repId=${repid}" \
-d "root='${root}'" \
$WORKING_FOLDER/src/07_SLiM/01_ph.slim

#--------------------------------------------------------------------------------

rm scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# Say done
echo "done"