#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=6
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-38


### THESE LISTS NEED TO FOLLOW THE SAME ORDER
INPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list1.txt)  ### A list of your R1 fastq files
INPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list2.txt)  ### A list of your R2 fastq files
OUTPUT1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list3.txt) ### A list of R1 output names -- just capture the meaningful part of the fastq names including R1 (remove the fq.gz or fq suffix)
OUTPUT2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" list4.txt) ### A list of R2 output names -- just capture the meaningful part of the fastq names including R2 (remove the fq.gz or fq suffix)

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF FASTP IN YOUR MACHINE/SERVER  = fastp v.0.20.1
module load fastp

### Trimming --  for each sample pair of raw fastq reads or for each library, we produce a pair of trimmed output files. 
### Remember to keep R1 and R2 in the output names created in lists 3 and 4 above
### This below uses 4 cores

fastp -w  4 -i $INPUT1 -I $INPUT2 -o $OUTPUT1\_trimmed.fastq.gz -O $OUTPUT2\_trimmed.fastq.gz