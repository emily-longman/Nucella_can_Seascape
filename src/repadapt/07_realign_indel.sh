#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-38

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "add_RG" ]
then echo "Working add_RG folder exist"; echo "Let's move on."; date
else echo "Working add_RG folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/add_RG; date
fi
