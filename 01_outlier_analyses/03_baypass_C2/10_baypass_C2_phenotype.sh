#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=baypass_C2_phenotype

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

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass in C2 mode using an omega file. 
# This option gives a constrast file under the core model. 
# The contrast in this case will be between the populations that are known good drillers (SBR, VD) with those that are bad drillers (SH, FC).

# Load modules 
module load gcc/13.3.0-xp3epyt
baypass=/gpfs1/home/e/l/elongman/software/baypass_public/sources/g_baypass

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/outlier_analyses/baypass

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "C2_Pheno" ]
then echo "Working C2_Pheno folder exist"; echo "Let's move on."; date
else echo "Working C2_Pheno folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/outlier_analyses/baypass/C2_Pheno; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/outlier_analyses/baypass/C2_Pheno

# Run baypass in C2 mode with a contrast file
$baypass -npop 19 \
-gfile $WORKING_FOLDER/outlier_analyses/baypass/genobaypass \
-poolsizefile $WORKING_FOLDER/outlier_analyses/baypass/poolsize \
-omegafile $WORKING_FOLDER/outlier_analyses/baypass/omega/NC_baypass_mat_omega.out \
-contrastfile $WORKING_FOLDER/outlier_analyses/guide_files/NC.baypass.C2_guide.file_Pheno_drilling.txt \
-d0yij 4 \
-outprefix NC_baypass_C2_Pheno \
-nthreads 20

