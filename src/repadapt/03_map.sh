#!/bin/bash
#SBATCH --partition=week #general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=02-00:00:00 #30:00:00
#SBATCH --cpus-per-task=15 #10
#SBATCH --mem-per-cpu=60G #30G
#SBATCH --array=10,30 #1-38

#--------------------------------------------------------------------------------
# My additions:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "bwa_output" ]
then echo "Working bwa_output folder exist"; echo "Let's move on."; date
else echo "Working bwa_output folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/bwa_output; date
fi

# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_bwa=https://depot.galaxyproject.org/singularity/bwa:0.7.17--h5bf99c6_8
repadapt_samtools=https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_0

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/Trim_map.txt
# Extract sample names/files
i=`awk -F "\t" '{print $6}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo "Sample i:" ${i}

INPUT1=${i}_R1_trimmed.fastq.gz
INPUT2=${i}_R2_trimmed.fastq.gz

# Change directory
cd $WORKING_FOLDER/data/processed/repadapt/bwa_output

# Map read
apptainer run $repadapt_bwa bwa mem -t 15 $WORKING_FOLDER/data/processed/repadapt/genome/N.canaliculata_assembly.fasta.softmasked.fa \
$WORKING_FOLDER/data/processed/repadapt/fastp/$INPUT1 \
$WORKING_FOLDER/data/processed/repadapt/fastp/$INPUT2  \
> $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.sam

# Build bam file (and remove sam file)
apptainer run $repadapt_samtools samtools view -Sb -q 10 --threads 15 $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.sam > $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.bam
rm $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.sam

# Sort bam file (remove unsorted bam file)
apptainer run $repadapt_samtools samtools sort --threads 15 $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.bam > $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\_sorted.bam
rm $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\.bam

# Index sorted bam file
apptainer run $repadapt_samtools samtools index $WORKING_FOLDER/data/processed/repadapt/bwa_output/$i\_sorted.bam

#--------------------------------------------------------------------------------

#### All these lists below need to follow same order
#INPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt) ### List of trimmed R1 fastq reads (produced with script 01)
#INPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt) ### List of trimmed R2 fastq reads (produced with script 01)
#OUTPUT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list3.txt) ### List of output names -- extract the meaningful part of the name from the trimmed reads names. There is a single output file for each pair of reads (each sample or library)

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF SAMTOOLS AND BWA IN YOUR MACHINE/SERVER  = bwa-mem v.0.7.17-r1188   samtools v.1.16.1
#module load bwa samtools

#### Here we map the trimmed reads to the ref genome and we produce a sam, then we convert it to bam, we sort it and finally we index it
#### We end up with 1 bam file per sample after this. If you had multiple libraries per sample, you'd end up with 1 bam per library


#bwa mem -t 4 Betula_pendula_subsp._pendula.fa $INPUT1 $INPUT2  > ./bwa_output/$OUTPUT\.sam
#samtools view -Sb -q 10 ./bwa_output/$OUTPUT\.sam > ./bwa_output/$OUTPUT\.bam
#rm ./bwa_output/$OUTPUT\.sam
#samtools sort --threads  4 ./bwa_output/$OUTPUT\.bam > ./bwa_output/$OUTPUT\_sorted.bam
#rm ./bwa_output/$OUTPUT\.bam
#samtools index ./bwa_output/$OUTPUT\_sorted.bam