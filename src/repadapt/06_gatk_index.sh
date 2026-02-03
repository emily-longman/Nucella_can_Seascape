#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x_%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=20G

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# This is to create gatk index of ref genome needed for indel realignment

# Change directory
cd $WORKING_FOLDER/data/processed/repadapt/genome

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0
repadapt_picard=https://depot.galaxyproject.org/singularity/picard:2.26.3--hdfd78af_0

# Index with Samtools
apptainer run $repadapt_samtools samtools faidx N.canaliculata_assembly.fasta.softmasked.fa

# Index with Picard
apptainer run $repadapt_picard picard -Xmx20G CreateSequenceDictionary \
R=N.canaliculata_assembly.fasta.softmasked.fa O=N.canaliculata_assembly.fasta.softmasked.dict

#--------------------------------------------------------------------------------

### This is to create gatk index of ref genome needed for indel realignment
### Keep the output files of this command in the same dir where you keep the reference genome fasta


###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF SAMTOOLS IN YOUR MACHINE/SERVER  = samtools v.1.16.1
#module load StdEnv/2020 samtools

#samtools faidx Betula_pendula_subsp._pendula.fa


#module purge
###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF PICARD IN YOUR MACHINE/SERVER  = Picard Tools v.2.26.3
#module load picard java


#java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=Betula_pendula_subsp._pendula.fa O=Betula_pendula_subsp._pendula.dict