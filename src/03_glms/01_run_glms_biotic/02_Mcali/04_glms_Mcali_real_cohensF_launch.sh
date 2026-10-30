#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=glms_Mcali_cohensF

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Submit job array
#SBATCH --array=1-996%100

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run the accompanying 0r_glms_Mcali_real_cohensF.R script. 
# For efficiency, there are 5 R scripts separating GLMs with the real data and chunks of the permutations. Note: walltimes vary depending on if running real data versus permutations.
# This will run glms and look at the association of the allele frequencies of the outlier SNPs and the Marine biotic data.

# Load modules 
module load gcc/13.3.0
module load R/4.4.1

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Script folder.
SCRIPT_FOLDER=$WORKING_FOLDER/src/03_glms/01_run_glms_biotic/02_Mcali

#--------------------------------------------------------------------------------

# Guide file 
guide_file=$WORKING_FOLDER/guide_files/scaffold_names_guide_file_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 18919, 2; 996 partitions each with 19 scaffold names)
# Scaffold name         # Partition/array
# Backbone_10001              1
# Backbone_10003              1
# Backbone_10004              1
# Backbone_10005              1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Change directory
cd $WORKING_FOLDER/data/processed/glms

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $guide_file | awk '{print $1}' > scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# List scaffold names
cat scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/glms

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "glms_chunk_analysis_Mcali" ]
then echo "Working glms_chunk_analysis_Mcali folder exist"; echo "Let's move on."; date
else echo "Working glms_chunk_analysis_Mcali folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/glms/glms_chunk_analysis_Mcali; date
fi

# Move to working directory
cd $WORKING_FOLDER/data/processed/glms/glms_chunk_analysis_Mcali

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "cohensF" ]
then echo "Working cohensF folder exist"; echo "Let's move on."; date
else echo "Working cohensF folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/glms/glms_chunk_analysis_Mcali/cohensF; date
fi

#--------------------------------------------------------------------------------

# Run R script
Rscript $SCRIPT_FOLDER/04_glms_Mcali_real_cohensF.R "${SLURM_ARRAY_TASK_ID}"

#--------------------------------------------------------------------------------

# Remove temporary scaffold guide file
rm scaffold.names.${SLURM_ARRAY_TASK_ID}.txt

# Say done
echo "done"