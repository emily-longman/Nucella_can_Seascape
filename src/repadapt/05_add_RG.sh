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
if [ -d "remove_RG" ]
then echo "Working remove_RG folder exist"; echo "Let's move on."; date
else echo "Working remove_RG folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/remove_RG; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_picard=https://depot.galaxyproject.org/singularity/picard:2.26.3--hdfd78af_0

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Trim_map.txt
# Extract sample names/files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i}
JAVAMEM=20G

# Add read groups, we start with our deduplicated bam files and we get a deduplicated bam with read groups assigned per sample/library
apptainer run $repadapt_picard picard -Xmx20G AddOrReplaceReadGroups \
I=$WORKING_FOLDER/data/processed/repadapt/remove_dup/${i}_sorted_dedup.bam \
O=$WORKING_FOLDER/data/processed/repadapt/remove_dup/${i}_sorted_dedup_RG.bam \
RGID=${i} RGLB=${i}\_LB RGPL=ILLUMINA RGPU=unit1 RGSM=${i}

#--------------------------------------------------------------------------------

### Keep the lists below with the same order
#INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt) ### List of input deduplicated bam files
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt) ### List of output names (just remove .bam) from the inputs
#NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list3.txt)   ### Here you want to extract the sample name from the input name, which is used to set the read IDs. Can be the same as list2.txt


### Here we add read groups, we start with our deduplicated bam files and we get a deduplicated bam with read groups assigned per sample/library
###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF PICARD IN YOUR MACHINE/SERVER  = Picard Tools v.2.26.3
#module load picard java

#java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=$INPUT O=$OUTPUT\_RG.bam RGID=$NAME RGLB=$NAME\_LB RGPL=ILLUMINA RGPU=unit1 RGSM=$NAME