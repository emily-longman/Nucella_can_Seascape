#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=30:00:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=100G
#SBATCH --array=1-2 #1-38

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "remove_dup" ]
then echo "Working remove_dup folder exist"; echo "Let's move on."; date
else echo "Working remove_dup folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/remove_dup; date
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

# Remove duplicates. We feed it a bam, and we get a deduplicated bam per sample/library
export JAVA_TOOL_OPTIONS="-Xmx18G"
apptainer run $repadapt_picard picard MarkDuplicates \
INPUT=$WORKING_FOLDER/data/processed/repadapt/bwa_output/${i}_sorted.bam \
OUTPUT=$WORKING_FOLDER/data/processed/repadapt/remove_dup/${i}_sorted\_dedup.bam \
METRICS_FILE=$WORKING_FOLDER/data/processed/repadapt/remove_dup/${i}_sorted.bam\_DUP_metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true

#--------------------------------------------------------------------------------

### The lists below need to follow the same order
#INPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt) ### list of input bam files (output of script 03)
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt) ### list of output names --  I usually just remove the suffix (.bam) to extract the name from the input, and then the command below will add _dedup.bam 

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF PICARD IN YOUR MACHINE/SERVER  = Picard Tools v.2.26.3
#module load picard java


### Here we remove duplicates. We feed it a bam, and we get a deduplicated bam per sample/library

#java -jar $EBROOTPICARD/picard.jar MarkDuplicates INPUT=$INPUT OUTPUT=$OUTPUT\_dedup.bam METRICS_FILE=$INPUT\_DUP_metrics.txt VALIDATION_STRINGENCY=SILENT REMOVE_DUPLICATES=true