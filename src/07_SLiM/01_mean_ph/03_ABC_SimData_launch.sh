#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=pH_ABC_SimData_launch

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=600G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying R script. 

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

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
if [ -d "ph_ABC" ]
then echo "Working ph_ABC folder exist"; echo "Let's move on."; date
else echo "Working ph_ABC folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/SLiM/ph_ABC; date
fi

#--------------------------------------------------------------------------------

# Run R script
Rscript --vanilla $WORKING_FOLDER/src/07_SLiM/01_mean_ph/03_ABC_SimData.R

#--------------------------------------------------------------------------------

# Say done
echo "done"