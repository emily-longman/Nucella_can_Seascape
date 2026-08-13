#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=variant_calling_bcftools

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=30:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=8G 

# Request CPU
#SBATCH --cpus-per-task=2

# Submit job array
#SBATCH --array=1-631%40

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%A_%a.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will use bcftools to call variants.
# Generate VCF containing genotype likelihoods for the bam file list.

# Note:
# Since the N. canaliculata genome consistent of many scaffolds, I chunked the genome into 50 scaffolds partitions to utilize an array structure.
# For each slurm array ID, bcftools is running over the 50 scaffolds using a loop then the vcf files are combined together.  

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

# This is the location where the reference genome is stored.
REFERENCE_FOLDER=/netfiles/pespenilab_share/Nucella/processed/Base_Genome/Base_Genome_Oct2024/Crassostrea_softmask

# This is the path to the reference genome.
REFERENCE=$REFERENCE_FOLDER/N.canaliculata_assembly.fasta.softmasked.fa

#--------------------------------------------------------------------------------

## Import master partition file 
GUIDE_FILE=$WORKING_FOLDER/fastq_to_vcf/guide_files/Scaffold_names_array.txt

#Example: -- the headers are just for descriptive purposes. The actual file has no headers. (dimensions: 2, 18919; 631 partitions)
# Scaffold name       # Partition
# Backbone_10001            1
# Backbone_10003            1
# Backbone_10004            1
# ....

#--------------------------------------------------------------------------------

# Determine partition to process 

# Echo slurm array task ID
echo ${SLURM_ARRAY_TASK_ID}

# Change directory 
cd $WORKING_FOLDER/fastq_to_vcf

# Using the guide file, extract the scaffold names associated based on the Slurm array task ID for a given partition
awk '$2=='${SLURM_ARRAY_TASK_ID}'' $GUIDE_FILE | awk '{print $1}' > genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}.txt

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf_bcftools" ]
then echo "Working vcf_bcftools folder exist"; echo "Let's move on."; date
else echo "Working vcf_bcftools folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools; date
fi

cd $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools

if [ -d "partitions" ]
then echo "Working partitions folder exist"; echo "Let's move on."; date
else echo "Working partitions folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions; date
fi

cd $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions

if [ -d "genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}" ]
then echo "Working genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID} folder exist"; echo "Let's move on."; date
else echo "Working genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID} folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}; date
fi

#--------------------------------------------------------------------------------

# Prepare bamlist
# This is a file with the name and full path of all the bam files to be processed.

# Move to bams folder
cd $WORKING_FOLDER/fastq_to_vcf/RGSM_final_bams

# Create bamlist for all Nucella samples
ls -d "$PWD/"*.bam > $WORKING_FOLDER/fastq_to_vcf/guide_files/Nucella_pops_bam_names_bcftools.list

BAMLIST=$WORKING_FOLDER/fastq_to_vcf/guide_files/Nucella_pops_bam_names_bcftools.list

#--------------------------------------------------------------------------------

# Use bcftools to call variants
# For the scaffolds in a given partition, generate a vcf file for each scaffold then cat them together. 

# The 'mpileup' part generates genotype likelihoods at each genomic position with coverage. The 'call' part makes the actual calls for the variant sites. 

cd $WORKING_FOLDER/fastq_to_vcf

# Cat file of scaffold names and start while loop
cat $WORKING_FOLDER/fastq_to_vcf/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}.txt | \
while read REGION 
do echo ${REGION}

bcftools mpileup -d 500 -C 50 -Ou -f $REFERENCE -b $BAMLIST -r $REGION -q 30 -Q 20 -a DP,AD | bcftools call -mv -Oz  > $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}/Nucella.bcftools.${REGION}.vcf.gz

done 

### bcftools mpileup parameters:
# -d --max-depth: At a position, read maximally INT reads per input file
# -C --adjust-MQ: Coefficient for downgrading mapping quality for reads containing excessive mismatches.
# -O --output-type b|u|z|v: output uncompressed bcf (u) - need this format for "bcftools call"
# -f --fasta-ref: The faidx-indexed reference file in FASTA format.
# -b --bam-list: list of input alignment files, one file per line.
# -r --regions: Only generate mpileup output in given regions. Requires the alignment files to be indexed.
# -q --min-MQ: Minimum mapping quality for a base to be considered
# -Q --min_BQ: Minimum base quality for a base to be considered
# -a --annotate: Comma-separated list of FORMAT and INFO tags to output.

### bcftools call parameters:
# -m --multiallelic-caller: alternative model for multiallelic and rare-variant calling designed to overcome known limitations in -c calling mode 
# -v --variants-only: output variant sites only
# -Oz: --output type b|u|z|v:output compressed vcf (z)

# The read depth should be adjusted to about twice the average read depth as higher read depths usually indicate problematic regions which are often enriched for artefacts

#--------------------------------------------------------------------------------

# Combine the vcf files from all of the scaffolds in a given partition

bcftools concat \
$WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}/*.vcf.gz \
-Oz -o $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}.vcf.gz

#--------------------------------------------------------------------------------

# Change directory
cd $WORKING_FOLDER/fastq_to_vcf

# Housekeeping
rm genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}.txt
rm -R $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/partitions/genome.scaffold.names.bcftools.${SLURM_ARRAY_TASK_ID}

#--------------------------------------------------------------------------------

date
echo "done"

