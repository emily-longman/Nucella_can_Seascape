#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_GWAS_phenotypic

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=3-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=30G 

# Request CPU
#SBATCH --cpus-per-task=25

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will perform the first step in running baypass. 

# Load modules 
module load gcc/13.3.0-xp3epyt
#module load netcdf-fortran/4.6.1 # Use to compile baypass
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/phenotypic_data

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "Baypass_GWAS" ]
then echo "Working Baypass_GWAS folder exist"; echo "Let's move on."; date
else echo "Working Baypass_GWAS folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/phenotypic_data/Baypass_GWAS; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/data/processed/phenotypic_data/Baypass_GWAS

# Run baypass - this will generate the omega file which will be used in the subsequent scripts
$baypass -npop 11 \
-gfile $WORKING_FOLDER/data/processed/phenotypic_data/baypass_input_files/genobaypass \
-poolsizefile $WORKING_FOLDER/data/processed/phenotypic_data/baypass_input_files/poolsize \
-omegafile $WORKING_FOLDER/data/processed/phenotypic_data/omega/NC_pheno_mat_omega.out \
-efile $WORKING_FOLDER/data/raw/phenotypic_data/Pheno_drilling_baypass_efile.txt \
-d0yij 4 \
-outprefix NC_pheno \
-npilot 100 -nthreads 25
