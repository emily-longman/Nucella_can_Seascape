#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=05:00:00
#SBATCH --mem-per-cpu=20G
#SBATCH --array=1-19

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "realign_indel_fix_RG" ]
then echo "Working realign_indel_fix_RG folder exist"; echo "Let's move on."; date
else echo "Working realign_indel_fix_RG folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/realign_indel_fix_RG; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_picard=https://depot.galaxyproject.org/singularity/picard:2.26.3--hdfd78af_0

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/samples_file.txt
# Extract sample names/files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i}
JAVAMEM=20G

# Add read groups, we start with our deduplicated bam files and we get a deduplicated bam with read groups assigned per sample/library
apptainer run $repadapt_picard picard -Xmx20G AddOrReplaceReadGroups \
I=$WORKING_FOLDER/data/processed/repadapt/realign_indel/${i}_sorted_dedup_RG_lanes_merged_realigned.bam \
O=$WORKING_FOLDER/data/processed/repadapt/realign_indel_fix_RG/${i}_sorted_dedup_RG_lanes_merged_realigned.bam \
RGID=${i} RGLB=${i}\_LB RGPL=ILLUMINA RGPU=unit1 RGSM=${i}
