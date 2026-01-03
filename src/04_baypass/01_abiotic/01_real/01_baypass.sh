#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_abiotic

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss
#SBATCH --time=3-00:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Request CPU
#SBATCH --cpus-per-task=20

# Submit job array
#SBATCH --array=1-9

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

# Guide files - Note: need to make sure the latter file is in the same order as the samples in the gfile
var_names=$WORKING_FOLDER/guide_files/Bio-oracle_enviro_vars_names.txt
guide_file=$WORKING_FOLDER/guide_files/Baypass_abiotic.txt

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/data/processed/GEA/baypass

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Extract enviro var 
var=`sed "${SLURM_ARRAY_TASK_ID}q;d" $var_names`
echo $var

# Using the guide file, extract the bio-oracle data associated based on the Slurm array task ID 
awk -v var="$var" '$1 == var { print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20}' $guide_file > ${var}.data.txt

# State data
cat ${var}.data.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/baypass

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "abiotic" ]
then echo "Working abiotic folder exist"; echo "Let's move on."; date
else echo "Working abiotic folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/baypass/abiotic; date
fi

cd $WORKING_FOLDER/data/processed/GEA/baypass/abiotic

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "${var}" ]
then echo "Working ${var} folder exist"; echo "Let's move on."; date
else echo "Working ${var} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/baypass/abiotic/${var}; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/GEA/baypass/abiotic/${var}

# Run baypass in aux covaraiate mode to estimate Bayes Factors
$baypass -npop 19 \
-gfile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/genobaypass \
-poolsizefile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/poolsize \
-omegafile $WORKING_FOLDER/data/processed/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out \
-efile $WORKING_FOLDER/data/processed/GEA/baypass/${var}.data.txt \
-d0yij 4 \
-auxmodel \
-outprefix NC_abiotic_${var} \
-nthreads 20

#--------------------------------------------------------------------------------

# Housekeeping
rm $WORKING_FOLDER/data/processed/GEA/baypass/${var}.data.txt

# Say done
echo "done"