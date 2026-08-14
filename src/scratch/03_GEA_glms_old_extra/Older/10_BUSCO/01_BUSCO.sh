#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=BUSCO

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=60G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script runs BUSCO on an assembly in the current directory. 
# As a result, you must cd to that directory before running.
#cd $WORKING_FOLDER/data/processed/BUSCO

# Also prior to script you must get the BUSCO lineage
#wget https://busco-data.ezlab.org/v5/data/lineages/mollusca_odb12.2025-07-01.tar.gz

#--------------------------------------------------------------------------------

# Call package
module load miniforge
#conda create --name busco_env_v6 #create and name the environment
conda activate busco_env_v6 #activate the environment
#conda install -c conda-forge -c bioconda busco=6.0.0 # install the program

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/data/processed/BUSCO

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "BUSCO" ]
then echo "Working BUSCO folder exist"; echo "Let's move on."; date
else echo "Working BUSCO folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/BUSCO; date
fi

#--------------------------------------------------------------------------------

# Reference
REFERENCE=$WORKING_FOLDER/data/raw/genome/N.canaliculata_assembly.fasta

# Run BUSCO for drosophila_odb12 lineage
LINEAGE=$WORKING_FOLDER/data/processed/BUSCO/lineage/mollusca_odb12

busco -m genome -i $REFERENCE -o N_canaliculata -l $LINEAGE

#--------------------------------------------------------------------------------

conda deactivate