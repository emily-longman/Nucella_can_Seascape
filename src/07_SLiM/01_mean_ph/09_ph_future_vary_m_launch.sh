#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ph_future_vary_m

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=01:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=2G

# Submit job array
#SBATCH --array=1-3

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

# Guide file 
GUIDE_FILE=$WORKING_FOLDER/guide_files/slim_ph_vary_migration.txt

# Determine parameters
m=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "m:"${m}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "ph_vary_m" ]
then echo "Working ph_vary_m folder exist"; echo "Let's move on."; date
else echo "Working ph_vary_m folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_vary_m; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_vary_m

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "m_${m}" ]
then echo "Working m_${m} folder exist"; echo "Let's move on."; date
else echo "Working m_${m} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_vary_m/m_${m}; date
fi

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_vary_m/m_${m}

# Set root
ROOT=$WORKING_FOLDER/data/processed/SLiM/ph_future/ph_vary_m/m_${m}

# Loop through iterations
for i in {1..500}
do

# Run slim script
slim \
    -d "repId=${i}" \
    -d "m=${m}" \
	-d "root='${ROOT}'" \
     $WORKING_FOLDER/src/07_SLiM/01_mean_ph/09_ph_future_vary_m.slim

done

#--------------------------------------------------------------------------------

# Say done
echo "done"