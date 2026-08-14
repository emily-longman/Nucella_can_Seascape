#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=haplotypecaller

# Specify partition
#SBATCH --partition=week

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=4-00:00:00 

# Request CPU
#SBATCH --cpus-per-task=10

# Submit job array
#SBATCH --array=1-19

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/HaplotypeCaller.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will call haplotypes (gVCFs) with GATK. 

# Load modules 
module load singularity
module load gatk/4.6.1.0
PICARD=/netfiles/nunezlab/Shared_Resources/Software/picard/build/libs/picard.jar
TABIX=/netfiles/nunezlab/Shared_Resources/Software/htslib/tabix
BGZIP=/netfiles/nunezlab/Shared_Resources/Software/htslib/bgzip

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

# Name of pipeline
PIPELINE=HaplotypeCaller

#--------------------------------------------------------------------------------

# Define parameters
CPU=10
echo "using #CPUs ==" $CPU
JAVAMEM=18G # Java memory

#--------------------------------------------------------------------------------

#Read Information
Group_library="Longman_2025"
Library_Platform="NovaseqX"
Group_platform="Longman_2025"

#--------------------------------------------------------------------------------

# Read guide files
# This is a file with the names of all the populations. 
GUIDE_FILE=$WORKING_FOLDER/guide_files/Populations.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Population
##   ARA
##   BMR
##   CBL
##   ...
##   VD

#--------------------------------------------------------------------------------

# Determine sample to process, "i" and read files
i=`awk -F "\t" '{print $1}' $GUIDE_FILE | sed "${SLURM_ARRAY_TASK_ID}q;d"`
echo $i

#--------------------------------------------------------------------------------

# This part of the pipeline will generate log files to record warnings and completion status

# Move to logs directory
cd $WORKING_FOLDER/logs

echo $PIPELINE

if [[ -e "${PIPELINE}.completion.log" ]]
then echo "Completion log exist"; echo "Let's move on."; date
else echo "Completion log doesnt exist. Let's fix that."; touch $WORKING_FOLDER/logs/${PIPELINE}.completion.log; date
fi

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "RGSM_final_bams" ]
then echo "Working RGSM_final_bams folder exist"; echo "Let's move on."; date
else echo "Working RGSM_final_bams folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/RGSM_final_bams; date
fi

if [ -d "haplotype_calling" ]
then echo "Working haplotype_calling folder exist"; echo "Let's move on."; date
else echo "Working haplotype_calling folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/haplotype_calling; date
fi

#--------------------------------------------------------------------------------

# Create a uniform read group with picard
java -Xmx$JAVAMEM -jar $PICARD AddOrReplaceReadGroups \
I=$WORKING_FOLDER/bams_merged/${i}.lanes_merged.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
RGLB=$Group_library \
RGPL=$Library_Platform \
RGPU=$Group_platform \
RGSM=${i}

#--------------------------------------------------------------------------------

# Index bam files
java -Xmx$JAVAMEM -jar $PICARD BuildBamIndex \
I=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam \
O=$WORKING_FOLDER/RGSM_final_bams/${i}.RG.bai

#--------------------------------------------------------------------------------

# Haplotype calling

singularity exec $GATK_SIF gatk --java-options "-Xmx${JAVAMEM}" HaplotypeCaller \
-R $REFERENCE \
-I $WORKING_FOLDER/RGSM_final_bams/${i}.RG.bam  \
-O $WORKING_FOLDER/haplotype_calling/${i}.g.vcf \
-ERC GVCF

#--------------------------------------------------------------------------------

# Compress and index with Tabix
$bgzip $WORKING_FOLDER/haplotype_calling/${i}.g.vcf
$tabix $WORKING_FOLDER/haplotype_calling/${i}.g.vcf.gz

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)