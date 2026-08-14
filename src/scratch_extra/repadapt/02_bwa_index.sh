#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5G

### Here we just create the bwa index of the reference genome, needed for mapping
### Keep the output of this command in the same dir where you keep the reference genome fasta

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_bwa=https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8

# Change directory
cd $WORKING_FOLDER/data/processed/repadapt/genome

apptainer run $repadapt_bwa bwa index -a bwtsw N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF BWA IN YOUR MACHINE/SERVER  = bwa-mem v.0.7.17-r1188
#module load bwa 


#bwa index -a bwtsw Betula_pendula_subsp._pendula.fa