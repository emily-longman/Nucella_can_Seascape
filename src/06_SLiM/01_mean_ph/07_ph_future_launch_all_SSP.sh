#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ph_future_all_SSP

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=02:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G

# Submit job array
#SBATCH --array=4 #1-5

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

# SSP scenario
GUIDE_FILE=$WORKING_FOLDER/guide_files/IPCC_pH_Scenarios.txt

# Extract one scenario
scenario=$(sed "${SLURM_ARRAY_TASK_ID}q;d" $GUIDE_FILE)
echo "IPCC Scenario: " ${scenario} 

#--------------------------------------------------------------------------------

# Generate Folders and files

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "ph_results_${scenario}" ]
then echo "Working ph_results_${scenario} folder exist"; echo "Let's move on."; date
else echo "Working ph_results_${scenario} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_results_${scenario}; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "per_simcycle" ]
then echo "Working per_simcycle folder exist"; echo "Let's move on."; date
else echo "Working per_simcycle folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_results_${scenario}/per_simcycle; date
fi

# Determine parameters

# Set root
ROOT=$WORKING_FOLDER/data/processed/SLiM/ph_future/ph_results_${scenario}/per_simcycle

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_future/ph_results_${scenario}/per_simcycle

# Loop through iterations
for i in {1..500}
do

# Run slim script
slim \
    -d "repId=${i}" \
	-d "root='${ROOT}'" \
     $WORKING_FOLDER/src/06_SLiM/01_mean_ph/07_ph_future_SSP${scenario}.slim

done

#--------------------------------------------------------------------------------

# Say done
echo "done"