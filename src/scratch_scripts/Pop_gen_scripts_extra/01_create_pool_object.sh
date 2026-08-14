#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=genotype_table

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=10:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will filter the vcf file created in the previous script. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# VCF is the path to the vcf file created in step 11 (note, must unzip file)
VCF=$WORKING_FOLDER/fastq_to_bam/vcf_freebayes/N_can_pops.vcf

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "pop_structure" ]
then echo "Working pop_structure folder exist"; echo "Let's move on."; date
else echo "Working pop_structure folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_structure; date
fi

cd pop_structure

if [ -d "genotype_table" ]
then echo "Working genotype_table folder exist"; echo "Let's move on."; date
else echo "Working genotype_table folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/pop_structure/genotype_table; date
fi

#--------------------------------------------------------------------------------

# Output genotype table (variable loci)
bcftools query -f '%CHROM %POS[\t%GT]\n' $VCF > $WORKING_FOLDER/pop_structure/genotype_table/N.can.pop.genotypes.txt

# -f : defines output format
# %CHROM : Print chromosome for each VCF line
# %POS: Print position for each VCF line
# [\t%GT]: FORMAT tags can be extracted using the square bracket [] operator. Prints the GT field
# \n : newline character