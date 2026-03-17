#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=300G
#SBATCH --array=1-38

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "realign_indel" ]
then echo "Working realign_indel folder exist"; echo "Let's move on."; date
else echo "Working realign_indel folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/realign_indel; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0
repadapt_gatk3=https://depot.galaxyproject.org/singularity/gatk:3.8--9

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Trim_map.txt
# Extract sample names/files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i}
JAVAMEM=20G

# Index with samtools
apptainer run $repadapt_samtools samtools index $WORKING_FOLDER/data/processed/repadapt/merge_lanes/${i}_sorted_dedup_RG.bam

# Generate a indel realigned bam file for each sample/library
apptainer run $repadapt_gatk3 gatk3 -Xmx20G -T RealignerTargetCreator \
-R $WORKING_FOLDER/data/processed/repadapt/genome/N.canaliculata_assembly.fasta.softmasked.fa \
-I $WORKING_FOLDER/data/processed/repadapt/merge_lanes/${i}_sorted_dedup_RG.bam \
-o $WORKING_FOLDER/data/processed/repadapt/realign_indel/${i}_sorted_dedup_RG\.intervals

echo "Finished Realigner Target Creator; Now, Indel Realigner"

apptainer run $repadapt_gatk3 gatk3 -Xmx20G -T IndelRealigner \
-R $WORKING_FOLDER/data/processed/repadapt/genome/N.canaliculata_assembly.fasta.softmasked.fa \
-I $WORKING_FOLDER/data/processed/repadapt/merge_lanes/${i}_sorted_dedup_RG.bam \
-targetIntervals $WORKING_FOLDER/data/processed/repadapt/realign_indel/${i}_sorted_dedup_RG\.intervals \
--consensusDeterminationModel USE_READS -o $WORKING_FOLDER/data/processed/repadapt/realign_indel/${i}_sorted_dedup_RG\_realigned.bam

#--------------------------------------------------------------------------------


### Keep the lists below with the same order
#INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt) ### list of bam files with read groups (output of script 05)
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)  ### list of output names (just remove .bam suffix from input list)

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF SAMTOOLS IN YOUR MACHINE/SERVER  = samtools v.1.16.1
#module load StdEnv/2020 samtools

#samtools index $INPUT


#module purge
#module load nixpkgs/16.09
###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF GATK IN YOUR MACHINE/SERVER  = GATK v.3.8
#module load java gatk/3.8

#### Here we generate a indel realigned bam file for each sample/library

#java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R Betula_pendula_subsp._pendula.fa -I $INPUT -o $OUTPUT\.intervals

#java -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R Betula_pendula_subsp._pendula.fa -I $INPUT -targetIntervals $OUTPUT\.intervals --consensusDeterminationModel USE_READS  -o $OUTPUT\_realigned.bam