#!/usr/bin/env bash

# In the command line, run the following command: sbatch path/to/this/file.sh

# Request cluster resources ----------------------------------------------------

# Name this job
#SBATCH --job-name=filter_vcf_LD

# Specify partition
#SBATCH --partition=general

# Request nodes
#SBATCH --nodes=1 

# Reserve walltime -- hh:mm:ss --30 hrs max
#SBATCH --time=2:00:00 

# Request memory for the entire job -- you can request --mem OR --mem-per-cpu
#SBATCH --mem=40G 

# Request CPU
#SBATCH --cpus-per-task=3

# Name output of this job using %x=job-name and %j=job-id
#SBATCH --output=./slurmOutput/%x.%j.out # Standard output

# Receive emails when job begins and ends or fails
#SBATCH --mail-type=ALL
#SBATCH --mail-user=emily.longman@uvm.edu 

#--------------------------------------------------------------------------------

# NOTE: LD filtering is impossible with pool-seq data! Pool-seq looses all LD information during sequencing.
# Also, plink is not compatible with pool data.

# This script will filter SNPs that are in LD. 

# Load software
plink=/gpfs1/home/e/l/elongman/software/plink_linux_x86_64_20241022/plink

#--------------------------------------------------------------------------------

# Define important file locations

# WORKING_FOLDER is the core folder where this pipeline is being run.
WORKING_FOLDER=/gpfs2/scratch/elongman/Nucella_can_Pop_Genomics/data/processed

#--------------------------------------------------------------------------------

# Generate Folders and files

# Move to working directory
cd $WORKING_FOLDER/fastq_to_vcf

# This part of the script will check and generate, if necessary, all of the output folders used in the script

if [ -d "vcf_clean_LD" ]
then echo "Working vcf_clean_LD folder exist"; echo "Let's move on."; date
else echo "Working vcf_clean_LD folder doesnt exist. Let's fix that."; mkdir $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD; date
fi

#--------------------------------------------------------------------------------

# Change directory 
cd $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD

# Create plink files (plink.bed, plink.bim and plink.fam) from VCF 
$plink --vcf $WORKING_FOLDER/fastq_to_vcf/vcf_clean/N.canaliculata_pops_filter.recode.vcf \
--allow-extra-chr --out $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink

# Reference the plink files 
$plink --bfile $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink \
--allow-extra-chr --no-sex --no-pheno --recode tab --out $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink

# Calculate r2 between all SNPs in dataset
$plink --file $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink \
--allow-extra-chr --r2 --out $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink_r2 --threads 3

# Make a list of SNPs with a r2 less than 0.8
$plink --file $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink \
--allow-no-sex --allow-extra-chr --indep-pairwise 100 10 0.8 --r2 \
--out $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink_indep_pairwise_100_10_0.8 --threads 3

# indep-pairwise requires 3 parameters:
# 1) a window size in variant count of kilobase
# 2) a variant count to shift the window at the end of each step
# 3) a pairwise r2 threshold - at each step, pairs of variants in the current window with squared correlation greater than the threshold are noted, 
# and variants are greedily pruned from the window until no such pairs remain

# Lastly, extract the pruned SNPs from the vcf file and create a new vcf that contains only SNPs not in LD
$plink --vcf $WORKING_FOLDER/fastq_to_vcf/vcf_freebayes/N.canaliculata_pops.vcf \
--recode vcf --allow-extra-chr \
--out $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink.LDfiltered_0.8 \
--extract $WORKING_FOLDER/fastq_to_vcf/vcf_clean_LD/N.canaliculata_pops_filter.recode.plink_indep_pairwise_100_10_0.8.prune.in