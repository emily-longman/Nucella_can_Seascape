#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_POD

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=3-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G 

# Request CPU
#SBATCH --cpus-per-task=20

# Submit job array
#SBATCH --array=1-5

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass on the psuedo-observed data set (PODs) to calibrate the XtX. 

# Load modules 
module load gcc/13.3.0-xp3epyt
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/outlier_analyses/baypass

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "POD" ]
then echo "Working POD folder exist"; echo "Let's move on."; date
else echo "Working POD folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/baypass/POD; date
fi

cd $WORKING_FOLDER/data/processed/outlier_analyses/baypass/POD

if [ -d "POD_${SLURM_ARRAY_TASK_ID}" ]
then echo "Working POD folder exist"; echo "Let's move on."; date
else echo "Working POD folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/outlier_analyses/baypass/POD/POD_${SLURM_ARRAY_TASK_ID}; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/outlier_analyses/baypass/POD/POD_${SLURM_ARRAY_TASK_ID}

# Run baypass on POD
$baypass -npop 19 \
-gfile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/G.NC.baypass.${SLURM_ARRAY_TASK_ID}.sim \
-poolsizefile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/poolsize \
-omegafile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out \
-d0yij 4 \
-outprefix NC_baypass_POD_${SLURM_ARRAY_TASK_ID} \
-nthreads 20

