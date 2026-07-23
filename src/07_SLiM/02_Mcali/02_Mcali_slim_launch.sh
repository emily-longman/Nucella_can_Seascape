#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=Mcali_slim_v4

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=10G 

# Submit job array
#SBATCH --array=1-936

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
if [ -d "Mcali_results_v4" ]
then echo "Working Mcali_results_v4 folder exist"; echo "Let's move on."; date
else echo "Working Mcali_results_v4 folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/Mcali_results_v4; date
fi

#--------------------------------------------------------------------------------

# Guide file 
GUIDE_FILE=$WORKING_FOLDER/guide_files/slim_Mcali_guide_file_v4.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
# Threshold      # k         # mag      # m         #N
# 1.0           0.05            1        0.005     5000
# ...           
# 2.05          0.20            1        0.005     5000

#--------------------------------------------------------------------------------

# Determine parameters
thresh=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
k=`awk -F "\t" '{print $2}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
mag=`awk -F "\t" '{print $3}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
m=`awk -F "\t" '{print $4}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
N=`awk -F "\t" '{print $5}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Threshold:"${thresh} "k:"${k} "mag:"${mag} "m:"${m} "N:"${N}

# Set root
ROOT=$WORKING_FOLDER/data/processed/SLiM/Mcali_results_v4

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/data/processed/SLiM/Mcali_results_v4

# Loop through iterations
for i in {1..30}
#for i in {26..40}
do

# Run slim script
slim \
	-d "thresh=${thresh}" \
    -d "k=${k}" \
    -d "mag=${mag}" \
    -d "m=${m}" \
    -d "N=${N}" \
    -d "repId=${i}" \
	-d "root='${ROOT}'" \
     $WORKING_FOLDER/src/07_SLiM/02_Mcali/02_Mcali.slim

done

#--------------------------------------------------------------------------------

# Say done
echo "done"