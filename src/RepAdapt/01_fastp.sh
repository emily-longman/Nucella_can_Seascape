#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-38

#--------------------------------------------------------------------------------
# My additions:

# copy reads to gpfs3tmp

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape
RAW_READS=/gpfs3tmp/pi/jcnunez/All_shortreads


# Make folders
cd $WORKING_FOLDER/data/processed
if [ -d "repadapt" ]
then echo "Working repadapt folder exist"; echo "Let's move on."; date
else echo "Working repadapt folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt; date
fi
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "fastp" ]
then echo "Working fastp folder exist"; echo "Let's move on."; date
else echo "Working fastp folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/fastp; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt=https://depot.galaxyproject.org/singularity/fastp:0.20.1--h8b12597_0

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Trim_map.txt
# Extract sample names/files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
INPUT1=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
INPUT2=`awk -F "\t" '{print $2}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i} "Read 1:" ${INPUT1} "Read 2:" ${INPUT2}

# Test run
apptainer run $repadapt fastp -w  4 -i $RAW_READS/$INPUT1 -I $RAW_READS/$INPUT2 -o $i\_R1_trimmed.fastq.gz -O $i\_R2_trimmed.fastq.gz

#--------------------------------------------------------------------------------

### THESE LISTS NEED TO FOLLOW THE SAME ORDER
#INPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt)  ### A list of your R1 fastq files
#INPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)  ### A list of your R2 fastq files
#OUTPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list3.txt) ### A list of R1 output names -- just capture the meaningful part of the fastq names including R1 (remove the fq.gz or fq suffix)
#OUTPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list4.txt) ### A list of R2 output names -- just capture the meaningful part of the fastq names including R2 (remove the fq.gz or fq suffix)

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF FASTP IN YOUR MACHINE/SERVER  = fastp v.0.20.1
#module load fastp

### Trimming --  for each sample pair of raw fastq reads or for each library, we produce a pair of trimmed output files. 
### Remember to keep R1 and R2 in the output names created in lists 3 and 4 above
### This below uses 4 cores

#fastp -w  4 -i $INPUT1 -I $INPUT2 -o $OUTPUT1\_trimmed.fastq.gz -O $OUTPUT2\_trimmed.fastq.gz