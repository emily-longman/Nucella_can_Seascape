#!/bin/bash
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --output=./slurmOutput/%x.%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 
#SBATCH --time=30:00:00
#SBATCH --mem-per-cpu=60G
#SBATCH --array=1 #1-631%100

#--------------------------------------------------------------------------------
# My additions/modifications:

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Seascape

# Make folders
cd $WORKING_FOLDER/data/processed/repadapt
if [ -d "mpileup_gz" ]
then echo "Working mpileup_gz folder exist"; echo "Let's move on."; date
else echo "Working mpileup_gz folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/data/processed/repadapt/mpilmpileup_gzeup; date
fi

###### CHANGE THE LINE OF CODE BELOW TO LOAD THE CORRECT VERSION OF BCFTOOLS IN YOUR MACHINE/SERVER  = bcftools v. 1.16
# Use apptainer to use the correct versions of programs
module load apptainer/1.3.4
repadapt_bcftools=https://depot.galaxyproject.org/singularity/bcftools:1.16--hfe4b78e_1

### If you have a very fragmented reference with 1000s of scaffolds, rather than sending 1000s of jobs, you can feed a list of chromosomes to the command rather than a single chromosome/scaffold name.
### So you can for example split 1000 scaffolds names into 10 lists of 100 scaffolds each, and feed those 10 lists to the command using a job array of size 10 (1-10) -- one list of scaffolds per job.

# Guide File
GUIDE_FILE=$WORKING_FOLDER/data/processed/repadapt/guide_files/chromosomes_guide_file.txt
# Determine partition to process 
# Echo slurm array task ID
echo "Slurm array:" ${SLURM_ARRAY_TASK_ID}
# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $GUIDE_FILE | awk '{print $1}' > chr.names.${SLURM_ARRAY_TASK_ID}.txt


REFERENCE=$WORKING_FOLDER/data/processed/repadapt/genome/N.canaliculata_assembly.fasta.softmasked.fa

#CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chromosomes.txt) ### list of chromosomes (can be found in the FASTA index file of reference genome -- .fai file). This species has 14, that's why array number is 14. We parallelize by chromosome. 

### Here we call SNPs. 

### list.txt is a list of ALL the realigned bam files of ALL samples.
# Note: made in 08.5_create_guide_file_part1

### If you have multiple bam per samples due to multiple libraries, merge them into one bam per sample before this step. - did this after step 5

### ploidymap.txt is a list that keeps the same order of samples in list.txt and should use the sample ID names you set in script 05. 
### For each sample ID name has the ploidy -- for a diploid species just create a tab separated file with 2 columns (first column: sample ID name set in script 05, second column: 2)

# Cat file of scaffold names and start while loop
cat chr.names.${SLURM_ARRAY_TASK_ID}.txt | \
while read scaffold 
do echo ${scaffold}

apptainer run $repadapt_bcftools bcftools mpileup -Ou -f $REFERENCE --bam-list samples.list.txt -q 5 -r $scaffold -I -a FMT/AD,FMT/DP | apptainer run $repadapt_bcftools bcftools call -S $WORKING_FOLDER/data/processed/repadapt/guide_files/ploidymap.txt -G - -f GQ -mv -Ov | gzip -c > $WORKING_FOLDER/data/processed/repadapt/mpileup_gz/$scaffold\.vcf.gz

done

# Housekeeping - remove intermediate files
rm chr.names.${SLURM_ARRAY_TASK_ID}.txt