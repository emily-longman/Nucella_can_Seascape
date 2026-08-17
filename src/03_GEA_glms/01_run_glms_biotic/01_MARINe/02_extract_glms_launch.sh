#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=extract_glms_biotic

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=100G 

# Submit job array
#SBATCH --array=1-14

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 03_extract_glms.R script. 

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# TMP_FOLDER - path to temporary folder with extra storage
TMP_FOLDER=/gpfs3tmp/pi/jcnunez

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/MARINe_vars_names.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. 
# Env variable name   
# B.gland_m
# B.gland_hm
# C.dalli_m
# C.dalli_hm
# M.cali_m
# M.cali_hm
# M.tross_m
# M.tross_hm
# N.can_m
# N.can_hm
# N.ost_m
# N.ost_hm
# P.och_m
# P.och_hm     

#--------------------------------------------------------------------------------

# Determine the environmental variable to process

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
i=`sed -n ${SLURM_ARRAY_TASK_ID}p $guide_file`

# State environmental variable
echo ${i}

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/GEA/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_per_env_var" ]
then echo "Working glms_per_env_var folder exist"; echo "Let's move on."; date
else echo "Working glms_per_env_var folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/GEA/glms/glms_per_env_var; date
fi

# Move to working directory
cd $TMP_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_per_env_var" ]
then echo "Working glms_per_env_var folder exist"; echo "Let's move on."; date
else echo "Working glms_per_env_var folder doesnt exist. Let's fix that."; mkdir $TMP_FOLDER/glms_per_env_var; date
fi

# Change directory
cd $TMP_FOLDER/glms_per_env_var

if [ -d "glms_${i}" ]
then echo "Working glms_${i} folder exist"; echo "Let's move on."; date
else echo "Working glms_${i} folder doesnt exist. Let's fix that."; mkdir $TMP_FOLDER/glms_per_env_var/glms_${i}; date
fi

#--------------------------------------------------------------------------------

# Run R script
Rscript --vanilla $WORKING_FOLDER/src/03_GEA_glms/01_run_glms_biotic/01_MARINe/02_extract_glms.R "${i}"

#--------------------------------------------------------------------------------

# Say done
echo "done"