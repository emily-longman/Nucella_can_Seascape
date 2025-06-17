#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=filter_vcf_bcftools

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=5:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=5G 

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# This script will filter the vcf file created in the previous script. 

# Load modules 
module load gcc/13.3.0-xp3epyt
module load bcftools/1.19-iq5mwek
module load vcftools/0.1.16

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf_bcftools_clean" ]
then echo "Working vcf_bcftools_clean folder exist"; echo "Let's move on."; date
else echo "Working vcf_bcftools_clean folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean; date
fi

#--------------------------------------------------------------------------------

# Filter vcf file

# Filter based on missing data (keep only variants in 75% of pops)
vcftools --gzvcf $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz \
--minQ 40 --maf 0.05 --max-missing 0.75 --recode --recode-INFO-all \
--out $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered

# --minQ: Includes only sites with Quality value above this threshold.
# --maf: Include only sites with a minor allele frequency great than or equal to this value.
# --max-missing: Exclude sites on the basis of the proportion of missing data (defined btween 0 and 1). 0 allows sites that are completely missing and 1 indicates no missing data.
# --recode --recode-INFO-all: Write a new VCF file that keeps all the INFO flags from the old vcf file.

# Generate a summary of the number of SNPs for each filter category (output suffix: .FILTER.summary)
vcftools --vcf $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered.recode.vcf \
--FILTER-summary --out $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered.recode.vcf
# Generate a file containing the depth per site summed across all individuals (output suffix: .idepth)
vcftools --vcf $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered.recode.vcf \
--depth --out $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered.recode.vcf


#--------------------------------------------------------------------------------

#### Filter with bcftools

# Check filtering parameters
#bcftools query -f 'QUAL>=30' $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz | sort | head -n 20

# Check number of high quality variants 
#bcftools filter -i 'QUAL>=30' $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz | grep -v -c '^#'
# Total number of variants 
#grep -v -c '^#' $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz 

# Filter high quality SNPs
#bcftools filter -i 'QUAL>=30' $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz \
#-Oz -o $WORKING_FOLDER/fastq_to_vcf/vcf_bcftools_clean/N.canaliculata_bcftools_pops_filtered.vcf.gz

# bcftools parameters:
# bcftools can filter-in or filter-out using the options -i and -e respectively
# -Oz --output type b|u|z|v: output compressed vcf (z)
# -o --output: write output to file with name as specified
# QUAL: Call quality 
# DP: Depth

### Example from https://samtools.github.io/bcftools/howtos/variant-calling.html
#bcftools filter -e'QUAL<10 || (RPB<0.1 && QUAL<15) || (AC<2 && QUAL<15) || MAX(DV)<=3 || MAX(DV)/MAX(DP)<=0.3' \
#$WORKING_FOLDER/fastq_to_vcf/vcf_bcftools/N.canaliculata_bcftools_pops.vcf.gz


#--------------------------------------------------------------------------------