#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=treemix

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=7-00:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will run baypass in core model (i.e., xtx mode) using an omega file. 

# Load modules 
# Call package (installed with conda)
module load python3.11-anaconda/2024.02-1
source ${ANACONDA_ROOT}/etc/profile.d/conda.sh
#conda create --name treemix #If you haven't already done so, create and name the environment
conda activate treemix #activate the environment
#conda install bioconda::treemix # If you haven't already done so, install the program

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# List of population names
POPS=$WORKING_FOLDER/pop_structure/guide_files/Treemix_pops.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/pop_structure

# This part of the script will check and generate, if necessary, all of the output folders used in the script
if [ -d "treemix" ]
then echo "Working treemix folder exist"; echo "Let's move on."; date
else echo "Working treemix folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_structure/treemix; date
fi

# Move to working directory
cd $WORKING_FOLDER/pop_structure/treemix

# Create output directory
if [ -d "treemix_output" ]
then echo "Working treemix_output folder exist"; echo "Let's move on."; date
else echo "Working treemix_output folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_structure/treemix/treemix_output; date
fi


#--------------------------------------------------------------------------------

# Reformat genobaypass file for Treemix input file 

# Change every space to a comma, then change every second comma to a space
sed -e "s/ /,/g" $WORKING_FOLDER/outlier_analyses/baypass/genobaypass | sed -E 's/(,[^,]*),/\1 /g' > $WORKING_FOLDER/pop_structure/treemix/Treemix.input

# Add population names to treemix input file
cat $POPS $WORKING_FOLDER/pop_structure/treemix/Treemix.input | gzip > $WORKING_FOLDER/pop_structure/treemix/Treemix.input.pop.names.gz

#--------------------------------------------------------------------------------

# Run treemix for 1 thru 19 migration events (1 per pop) with blocks of 1000 SNPs

# Do a for loop to represent 19 migration events; use 1000 snp blocks
for i in {0..19};
do  treemix -i Treemix.input.pop.names.gz -se -k 1000 -m $i -o $WORKING_FOLDER/pop_structure/treemix/treemix_output/Treemix.Output.$i ;
done

# -se: calculate standard errors of migration weights
# -k: number of SNPs per block for estimation of covariance matrix; accounts for the fact that nearby SNPs are not independent
# -m: number of migration edges

#--------------------------------------------------------------------------------

# Deactivate conda
conda deactivate