#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=ph_slim_expandedparam

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Submit job array
#SBATCH --array=801-960 #1-800 #1-960

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
if [ -d "ph_results_expandedparam" ]
then echo "Working ph_results_expandedparam folder exist"; echo "Let's move on."; date
else echo "Working ph_results_expandedparam folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_results_expandedparam; date
fi

#--------------------------------------------------------------------------------

# Guide file 
GUIDE_FILE=$WORKING_FOLDER/guide_files/slim_ph_guide_file_expandedparam.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
# Threshold      # K_1         # K_2     # m
# 7.7             0              0       0.01  
# ...           
# 8.4           0.001           0.001    0.0001

#--------------------------------------------------------------------------------

# Determine parameters
thresh=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
k_1=`awk -F "\t" '{print $2}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
k_2=`awk -F "\t" '{print $3}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
m=`awk -F "\t" '{print $4}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Threshold:"${thresh} "k_1:"${k_1} "k_2:"${k_2} "m:"${m}

# Set root
ROOT=$WORKING_FOLDER/data/processed/SLiM/ph_results_expandedparam

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/ph_results_expandedparam

# Loop through iterations (6-30 for 801-960 so do 20-30 for 1-800)
for i in {6..30}
do

# Run slim script
slim \
	-d "thresh=${thresh}" \
    -d "k_1=${k_1}" \
    -d "k_2=${k_2}" \
    -d "m=${m}" \
    -d "repId=${i}" \
	-d "root='${ROOT}'" \
     $WORKING_FOLDER/src/07_SLiM/02_ph_expandedparam.slim

done

#--------------------------------------------------------------------------------

# Say done
echo "done"