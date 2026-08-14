#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_xtx_genomic_offset

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=3-00:00:00

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
echo "Doing Baypass xtx run:" ${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/baypass

# Generate folder
if [ -d "xtx" ]
then echo "Working xtx folder exist"; echo "Let's move on."; date
else echo "Working xtx folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/baypass/xtx; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/baypass/xtx

# Run baypass in standard covariate mode to estimate Bayes Factors
$baypass -npop 19 \
-gfile $WORKING_FOLDER/data/processed/baypass/input_files/genobaypass \
-poolsizefile $WORKING_FOLDER/data/processed/baypass/input_files/poolsize \
-omegafile $WORKING_FOLDER/data/processed/baypass/omega/NC_baypass_mat_omega.out \
-d0yij 4 \
-outprefix NC_run${SLURM_ARRAY_TASK_ID} \
-nthreads 20

#--------------------------------------------------------------------------------

# Say done
echo "done"