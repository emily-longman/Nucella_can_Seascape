#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_Mtross_9pop

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=2-05:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=50G 

# Request CPU
#SBATCH --cpus-per-task=20

# Submit job array
#SBATCH --array=1-5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass for a given abiotic variable.

# Load modules 
module load gcc/13.3.0-xp3epyt
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Echo slurm array task ID
echo "Doing Baypass run:" ${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/baypass

# Generate folder
if [ -d "biotic" ]
then echo "Working biotic folder exist"; echo "Let's move on."; date
else echo "Working biotic folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/baypass/biotic; date
fi

# Change directory
cd $WORKING_FOLDER/data/processed/baypass/biotic

# Generate folder
if [ -d "Mtross_mean_9pop" ]
then echo "Working Mtross_mean_9pop folder exist"; echo "Let's move on."; date
else echo "Working Mtross_mean_9pop folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/baypass/biotic/Mtross_mean_9pop; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/baypass/biotic/Mtross_mean_9pop

# Run baypass in standard covariate mode to estimate Bayes Factors
$baypass -npop 9 \
-gfile $WORKING_FOLDER/data/processed/baypass/input_files/subset9pop.genobaypass \
-poolsizefile $WORKING_FOLDER/data/processed/baypass/input_files/subset9pop.poolsize \
-omegafile $WORKING_FOLDER/data/processed/baypass/omega_subset9pop/NC_subset_baypass_mat_omega.out \
-efile $WORKING_FOLDER/guide_files/Baypass_biotic_Mtross_9pop.txt \
-d0yij 4 \
-outprefix NC_biotic_Mtross_mean_run${SLURM_ARRAY_TASK_ID} \
-nthreads 20

#--------------------------------------------------------------------------------

# Say done
echo "done"