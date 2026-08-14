#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=build_DB_GATK

# Specify partition
#SBATCH --partition=bigmem

# Request nodes
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=20:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=150G 

# Submit job array
#SBATCH --array=1-38%20

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will merge gVCFs into a unified database for genotype calling. 

# Load modules 
module load singularity
module load gatk/4.6.1.0

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

# Name of pipeline
PIPELINE=Build_DB

#--------------------------------------------------------------------------------

# Define parameters
CPU=$SLURM_CPUS_ON_NODE
echo "using #CPUs ==" $SLURM_CPUS_ON_NODE
JAVAMEM=150G # Java memory

#--------------------------------------------------------------------------------

# Read guide files

# This is a file with the sample names and the path to the sample vcf per line. 
SAMPLE_MAP=$WORKING_FOLDER/guide_files/Samples_to_haplotype.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers.
## Population  g.vcf name
##   ARA      /gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam/haplotype_calling/ARA.g.vcf.gz 
##   BMR      /gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam/haplotype_calling/BMR.g.vcf.gz
##   CBL      /gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam/haplotype_calling/CBL.g.vcf.gz
##   ...
##   VD       /gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed/fastq_to_bam/haplotype_calling/VD.g.vcf.gz


# This is a file with the contig/scaffold names of the genome broken up into 38 partitions/chunks with 500 contigs/scaffolds listed. 
# We will utilize an array to work across the partitions.  
DB_INTERVALS=$WORKING_FOLDER/guide_files/DB_intervals_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 18919, 2 ; 38 partitions)
## Scaffold name         Partition
## Backbone_10001            1
## Backbone_10003            1
## ...
## ntLink_9                  38

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $DB_Intervals | awk '{print $1}' > $WORKING_FOLDER/guide_files/scaffold.names.${SLURM_ARRAY_TASK_ID}.list

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

if [ -d "DB_Nucella_canaliculata" ]
then echo "Working DB_Nucella_canaliculata folder exist"; echo "Let's move on."; date
else echo "Working DB_Nucella_canaliculata folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/DB_Nucella_canaliculata; date
fi

if [ -d "tmp_mergevcf.${SLURM_ARRAY_TASK_ID}" ]
then echo "Working tmp_mergevcf folder exist"; echo "Let's move on."; date
else echo "Working tmp_mergevcf folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/tmp_mergevcf; date
fi

#--------------------------------------------------------------------------------

# Merge VCFs using GenomicsDBImport
singularity exec $GATK_SIF gatk --java-options "-Xmx${JAVAMEM}" GenomicsDBImport \
--genomicsdb-workspace-path $WORKING_FOLDER/DB_Nucella_canaliculata \
--batch-size 50 \
--sample-name-map $SAMPLE_MAP \
--tmp-dir=$WORKING_FOLDER/tmp_mergevcf.${SLURM_ARRAY_TASK_ID} \
--reader-threads $CPU \
-L $WORKING_FOLDER/guide_files/scaffold.names.${SLURM_ARRAY_TASK_ID}.list

# --genomicsdb-workspace-path: workspace for GenomicsDB
# --batch-size: batch size controls the number of samples for which readers are open at once and therefore provides a way to minimize memory consumption
# --sample-name-map: is tab deliminated txt file; using it saves the tool from having to download the GVCF headers in order to determine the sample names
# --temp-dir: temporary directory to use
# -L: specifies the genomic intervals over which to operate

#--------------------------------------------------------------------------------

# Housekeeping

rm $WORKING_FOLDER/guide_files/scaffold.names.${SLURM_ARRAY_TASK_ID}.list
rm -r $WORKING_FOLDER/tmp_mergevcf.${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

# Inform that sample is done

# This part of the pipeline will notify the completion of run i. 

echo ${i} " completed" >> $WORKING_FOLDER/logs/${PIPELINE}.completion.log

echo "pipeline completed" $(date)